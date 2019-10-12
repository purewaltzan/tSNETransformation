# tSNETransformation
This library is a plugin of Seurat for t-SNE transformation. Users can do t-SNE transformation according to steps below:
1. Download library. All files should be in the same folder.

2. Running R script:<br>
```
library(seurat)
source('tsnetplugin.R') # t-SNE transformation library should be located here
mca <- CreatSeuratObject(counts = yourdata, project = yourProjectName)
mca <- NormalizeData(object = mca, normalization.method='RC')
mca <- tsneTransform(object = mca, 
                     perp = 30, # Perplexity in t-SNE 
                     dim = 30, # Dimension of space for t-SNE transformation after PCA 
                     platform = ['R'|'MATLAB'], # Select software to do t-SNE transformation. Suggest using MATLAB for matrix with large number of cells 
                     )
```

3. Log Normalization expression matrix for Heatmap and FeaturePlot <br>
```
mca <- NormalizeData(object = mca, normalization.method='LogNormalization',scale.factor = 10000) 
mca <- FindVariableFeatures(object = mca,...) 
mca <- ScaleData(object = mca,...) 
```

4. Further analysis steps in Seurat. such as 
```
mca <- FindNeibour(mca,...)
...
```
or
```
mca <- RunUMAP(mca, reduction = 'tsnetransform', dim=1:30,...) 
...
```

Important: Switch of step 2 and step 3 is legal. But t-SNE transformation do not keep the preference of independent genes.
