import os
import re
import time
import logging
import argparse
import datetime
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from transformers import T5Tokenizer, T5Model
import wandb
from tqdm import tqdm
import csv

from train_utils import load_model, process_sequence
from datasets import SequenceDataset


def train(model, tokenizer, dataloader, criterion, optimizer, device, epoch, run, csv_writer, csv_file):
    model.train()
    total_loss = 0
    num_batches = len(dataloader)
    logging.info(f'Starting training for epoch {epoch+1}. Total batches: {num_batches}')
    for batch_idx, batch in tqdm(enumerate(dataloader), desc=f'Epoch {epoch+1}', total=num_batches):
        start_time = time.time()

        seq1_list, seq2_list, labels = batch  # Unpack the batch

        # Process sequences to replace invalid amino acids
        seq1_list = [process_sequence(seq) for seq in seq1_list]
        seq2_list = [process_sequence(seq) for seq in seq2_list]

        # Tokenize sequences with padding and truncation
        inputs1 = tokenizer(seq1_list, return_tensors='pt', padding=True, truncation=True, max_length=128)
        inputs2 = tokenizer(seq2_list, return_tensors='pt', padding=True, truncation=True, max_length=128)

        # Move inputs and labels to the device
        inputs1 = {key: val.to(device) for key, val in inputs1.items()}
        inputs2 = {key: val.to(device) for key, val in inputs2.items()}
        labels = labels.to(device)

        # Convert labels from 0/1 to -1/1 for CosineEmbeddingLoss
        labels = labels * 2 - 1  # Now labels are -1 or 1

        # Forward pass through the encoder only
        outputs1 = model.encoder(input_ids=inputs1['input_ids'], attention_mask=inputs1['attention_mask'])
        outputs2 = model.encoder(input_ids=inputs2['input_ids'], attention_mask=inputs2['attention_mask'])

        # Obtain embeddings by mean pooling over the sequence length
        embeddings1 = outputs1.last_hidden_state.mean(dim=1)  # Shape: (batch_size, hidden_size)
        embeddings2 = outputs2.last_hidden_state.mean(dim=1)

        # Compute the cosine embedding loss
        loss = criterion(embeddings1, embeddings2, labels.float())

        # Backpropagation and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        loss = loss.item()
        total_loss += loss

        # Log the loss per batch to wandb
        run.log({'Train Loss per Batch': loss}, step=epoch * len(dataloader) + batch_idx)
        csv_writer.writerow(['train',epoch + 1, batch_idx + 1, loss])
        print('train',epoch + 1, batch_idx + 1, loss)
        csv_file.flush()  # Ensure the data is written to the file in real time

        end_time = time.time()
        batch_time = end_time - start_time

        if batch_idx % 10 == 0:
            logging.debug(f'Epoch [{epoch+1}], Batch [{batch_idx+1}/{num_batches}], Loss: {loss:.4f}, Time per batch: {batch_time:.2f}s')

    avg_loss = total_loss / len(dataloader)
    logging.info(f'Epoch [{epoch+1}] Training completed. Average Loss: {avg_loss:.4f}')
    return avg_loss

def validate(model, tokenizer, dataloader, criterion, device, epoch, run, csv_writer):
    model.eval()
    total_loss = 0
    num_batches = len(dataloader)
    logging.info(f'Starting validation for epoch {epoch+1}. Total batches: {num_batches}')
    with torch.no_grad():
        for batch_idx, batch in enumerate(dataloader):
            start_time = time.time()

            seq1_list, seq2_list, labels = batch  # Unpack the batch

            # Process sequences
            seq1_list = [process_sequence(seq) for seq in seq1_list]
            seq2_list = [process_sequence(seq) for seq in seq2_list]

            # Tokenize sequences with padding and truncation
            inputs1 = tokenizer(seq1_list, return_tensors='pt', padding=True, truncation=True)
            inputs2 = tokenizer(seq2_list, return_tensors='pt', padding=True, truncation=True)

            # Move inputs and labels to the device
            inputs1 = {key: val.to(device) for key, val in inputs1.items()}
            inputs2 = {key: val.to(device) for key, val in inputs2.items()}
            labels = labels.to(device)

            # Convert labels from 0/1 to -1/1 for CosineEmbeddingLoss
            labels = labels * 2 - 1  # Now labels are -1 or 1

            # Forward pass through the encoder only
            outputs1 = model.encoder(input_ids=inputs1['input_ids'], attention_mask=inputs1['attention_mask'])
            outputs2 = model.encoder(input_ids=inputs2['input_ids'], attention_mask=inputs2['attention_mask'])

            # Obtain embeddings by mean pooling over the sequence length
            embeddings1 = outputs1.last_hidden_state.mean(dim=1)
            embeddings2 = outputs2.last_hidden_state.mean(dim=1)

            # Compute the cosine embedding loss
            loss = criterion(embeddings1, embeddings2, labels.float())

            total_loss += loss.item()

            # Log the loss per batch to wandb
            run.log({'Validation Loss per Batch': loss.item()}, step=epoch * len(dataloader) + batch_idx)
            csv_writer.writerow(['validation',epoch + 1, batch_idx + 1, loss.item()])

            end_time = time.time()
            batch_time = end_time - start_time

            if batch_idx % 10 == 0:
                logging.debug(f'Epoch [{epoch+1}], Validation Batch [{batch_idx+1}/{num_batches}], Loss: {loss.item():.4f}, Time per batch: {batch_time:.2f}s')

    avg_loss = total_loss / len(dataloader)
    logging.info(f'Epoch [{epoch+1}] Validation completed. Average Loss: {avg_loss:.4f}')
    return avg_loss

def main():
    parser = argparse.ArgumentParser(description='Train protein sequence similarity model')
    parser.add_argument('--verbosity', type=int, default=1, help='Verbosity level for logging (0 = minimal, higher = more logs)')
    parser.add_argument('--batch_size', type=int, default=2, help='Batch size for training and validation')
    parser.add_argument('--epochs', type=int, default=1, help='Number of epochs to train')
    args = parser.parse_args()

    # Set up logging level based on verbosity
    logging_level = logging.INFO if args.verbosity == 1 else logging.DEBUG if args.verbosity > 1 else logging.WARNING
    logging.basicConfig(level=logging_level, format='%(asctime)s - %(levelname)s - %(message)s')

    # Save model and CSV in a timestamped directory
    model_save_dir = f'models/{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
    os.makedirs(model_save_dir, exist_ok=True)
    csv_file_path = os.path.join(model_save_dir, 'losses.csv')

    # Initialize the CSV file for logging
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Setting','Epoch', 'Average Train Loss', 'Average Validation Loss'])

        # Log start of main
        logging.info('Starting main function')

        # Initialize wandb in online mode to enable logging and visualization
        run = wandb.init(mode='offline', project='protein_sequence_similarity', name='Training Run')

        # Load model and tokenizer
        logging.info('Loading model and tokenizer')
        start_time = time.time()
        model, tokenizer = load_model()
        end_time = time.time()
        logging.info(f'Loaded model and tokenizer in {end_time - start_time:.2f}s')

        # Set device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model.to(device)
        logging.info(f'Using device: {device}')

        batch_size = args.batch_size
        logging.info(f'Batch size: {batch_size}')

        # Create training dataset and dataloader
        logging.info('Loading training dataset')
        train_dataset = SequenceDataset(file_path='/scratch/gpfs/vv7118/projects/vcmsa/data/ProtTucker/train66k.fasta', n=5000)
        train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        end_time = time.time()
        logging.info(f'Loaded training dataset in {end_time - start_time:.2f}s')

        # Create validation dataset and dataloader
        logging.info('Loading validation dataset')
        start_time = time.time()
        val_dataset = SequenceDataset(file_path='/scratch/gpfs/vv7118/projects/vcmsa/data/ProtTucker/val200.fasta')
        val_dataloader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        end_time = time.time()
        logging.info(f'Loaded validation dataset in {end_time - start_time:.2f}s')

        criterion = nn.CosineEmbeddingLoss()
        optimizer = optim.Adam(model.parameters(), lr=1e-5)

        num_epochs = args.epochs  # Set the number of epochs
        logging.info(f'Number of epochs: {num_epochs}')

        best_val_loss = float('inf')  # Initialize best validation loss

        for epoch in range(num_epochs):
            logging.info(f'Starting epoch {epoch+1}/{num_epochs}')

            # Training loop
            avg_train_loss = train(model, tokenizer, train_dataloader, criterion, optimizer, device, epoch, run, csv_writer, csv_file)

            # Validation loop
            avg_val_loss = validate(model, tokenizer, val_dataloader, criterion, device, epoch, run, csv_writer)

            logging.info(f"Epoch [{epoch+1}/{num_epochs}], Train Loss: {avg_train_loss:.4f}, Val Loss: {avg_val_loss:.4f}")

            # Log the average losses and epoch number to wandb
            run.log({'Epoch': epoch+1, 'Average Train Loss': avg_train_loss, 'Average Validation Loss': avg_val_loss})

            # Write losses to CSV
            # csv_writer.writerow([epoch + 1, avg_train_loss, avg_val_loss])

            # Check if the current validation loss is the best
            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                # Save the best model
                torch.save(model.state_dict(), os.path.join(model_save_dir, 'best_model.pt'))
                logging.info(f"Best model saved at epoch {epoch+1} with validation loss {avg_val_loss:.4f}")

        # Finish the wandb run
        run.finish()
        logging.info('Training completed')

if __name__ == '__main__':
    main()