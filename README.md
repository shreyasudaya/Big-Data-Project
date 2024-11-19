# Covid-19 variant classification using GNN
By Shreyas Udaya, Harsha N.P., Kiran Reddy R.

![image](../fig/vig.png)

## Requirements
Pytorch,
Seaborn,
Torch Geometric,
Matplotlib,
Sklearn,
Torchvision

## Dataset

The Dataset consists of heatmap images split into the image-net format i.e. split into the train and test directories in the ratio of 3:1 for train and test.
The Dataset can be found at https://www.kaggle.com/datasets/shreyasudaya/fcgr-covid-variants 

## Evaluation and Train Code

Data preparation follows the algorithm

- Run code:
    1. Take gnn-classifier.ipynb either on colab or kaggle
    2. Run the notebook and convert runtime type to GPU

## Directory Structure

Dataset is also available in directory Dataset, with it consisting of kmers, one-hot and fcgr feature extractions. Alongside this is Feature Ext, which contains the files in order to extract the feature and convert it.

```
project
│   README.md
│   gnn_classifier.ipynb
|   splitset.py    
│
└───Dataset
│   │   
│   └───FCGR
|   └───kmer  
│   └───one-hot    
│        
└───Feature Ext
    │   kprocheatmap.py
    │   kproc.py
    |   proc.py
```
