Ablation Study for vcMSA (by Wuwei)
=========

### Ablation on Embedding Extraction Layers

* Run `./layer_experiment/layer_experiment_all.sh` to get the original vcMSA results using concatenated amino acid embeddings from last 16 layers.
* Run `./layer_experiment/layer_experiment.sh` to get vcMSA results using amino acid embeddings from different model layers, one at a time.
* Run `./layer_experiment/eval.sh` to get performance evaluation results, using sum-of-pairs score and total column score as metrics.
* Run `./layer_experiment/visualize.ipynb` to get result visualizations.

### Ablation on Pooling Methods

* Run `./pooling_experiment/pooling_experiment.sh` to get vcMSA results using different pooling methods to obtain sequence embeddings from amino acid embeddings.
* Run `./pooling_experiment/eval.sh` to get evaluation results, using sum-of-pairs score and total column score as metrics.
* Run `./pooling_experiment/visualize.ipynb` to get result visualizations.
