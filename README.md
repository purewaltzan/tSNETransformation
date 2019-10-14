# tSNETransformation
This library is a plugin of Seurat for t-SNE transformation. Users can do t-SNE transformation according to steps below:
1. Download library. All files should be in the same folder.

2. Before do t-SNE transformation, ribosomal protein genes are strongly suggested being deleted. The code for mouse gene is:
```
featuress <- grep(pattern = 'Rps', x = rownames(yourdata.original), value = TRUE)
featuresl <- grep(pattern = 'Rpl', x = rownames(yourdata.original), value = TRUE)
yourdata <- yourdata.original[setdiff(setdiff(rownames(yourdata.original),featuress),featuresl),]
```

3. Running R script:<br>
```
library(seurat)
source('tsnetplugin.R') # set the path to the location of library you downloaded
mca <- CreatSeuratObject(counts = yourdata, project = yourProjectName)
mca <- NormalizeData(object = mca, normalization.method='RC')
mca <- tsneTransform(object = mca, 
                     perp = 30, # Perplexity in t-SNE 
                     dim = 30, # Dimension of space for t-SNE transformation after PCA 
                     platform = ['R'|'MATLAB'], # Select software to do t-SNE transformation. Suggest using MATLAB for matrix with large number of cells 
                     )
```

4. Log Normalization expression matrix for Heatmap and FeaturePlot <br>
```
mca <- NormalizeData(object = mca, normalization.method='LogNormalization',scale.factor = 10000) 
mca <- FindVariableFeatures(object = mca,...) 
mca <- ScaleData(object = mca,...) 
```

5. Further analysis steps in Seurat. such as 
```
mca <- FindNeibour(mca,reduction='tsnetransformation',...)  # Please annouce reduction='tsnetransformation', or the clustering will follow the results of PCA
...
```
or
```
mca <- RunUMAP(mca, reduction = 'tsnetransform', dim=1:30,...)  # Please annouce reduction='tsnetransformation', or the clustering will follow the results of PCA
...
```

Important: Switch of step 3 and step 4 is legal. But t-SNE transformation do not keep the preference of independent genes.
