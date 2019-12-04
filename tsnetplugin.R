if(!require('tsne')){
  install.packages('tsne')
  library(tsne)
}
if(!require('rstudioapi')){
  install.packages('rstudioapi')
  library(rstudioapi)
}
if(!require('matlabr')){
  install.packages('matlabr')
  library(matlabr)
}
libpath<-dirname(rstudioapi::getSourceEditorContext()$path)

dim.reduction <- setClass(
  Class = "DimReduc",
  slots = c(
    cell.embeddings = 'matrix',
    feature.loadings = 'matrix',
    feature.loadings.projected = 'matrix',
    assay.used = 'character',
    stdev = 'numeric',
    key = 'character',
    jackstraw = 'JackStrawData',
    misc = 'list'
  )
)

tsneTransform<-function(object,
                        perp=30,
                        dim=30,
                        platform='R',
                        reduction = "tsnetransform",
                        assay = NULL,
                        weight.by.var = TRUE,
                        verbose = TRUE,
                        reduction.key = "tsnet_",
                        seed.use = 42,
                        approx = TRUE,
                        outputfile = NULL,
                        deleteTemp = FALSE,
                        slot = 'data',
                        ...){
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  cells <- colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  raw_data<-as.matrix(GetAssayData(object,slot))
  if(length(object@assays$RNA@var.features>10)){
    raw_data<-raw_data[object@assays$RNA@var.features,]
  }
  if(platform=='outputfiles'){
    if(is.null(outputfile)){
      outputfile<-paste(libpath,'/tsnetoutput.csv',sep='')
      tsner<-read.csv(outputfile,header = FALSE, row.names = NULL)
    }
  }else{
    if(platform=='R'){
      warning("Trying tSNE transformation in R. If there are more than 2000 cells, it will take very long time. Try to use platform='MATLAB'")
      tsner<-tsne(t(raw_data),
                perplexity = perp,
                initial_dims = dim,
                k = dim,
                )
      }else{
        if(platform=='MATLAB'){
          get_matlab()
          if(have_matlab()){
            write.csv(t(raw_data),file=paste(libpath,'/data.csv',sep=''))
            write.csv(c(perp,dim),file=paste(libpath,'/para.csv',sep=''))
            run_matlab_script(fname=paste(libpath,'/tsnet.m',sep=''), verbose = TRUE, desktop = FALSE,
                          splash = FALSE, display = FALSE, wait = TRUE,
                          single_thread = FALSE, ...)
            outputfile<-paste(libpath,'/tsnetoutput.csv',sep='')
            tsner<-read.csv(outputfile,header = FALSE, row.names = NULL)
            #deleteTemp = TRUE
            }else{
              warning('There is no MATLAB detected in this environment, attempting to use output files')
              if(is.null(outputfile)){
                outputfile<-paste(libpath,'/tsnetoutput.csv',sep='')
                }
              tsner<-try(outputfile,header = FALSE, row.names = NULL)
              if(!is.array(tsner)){
                stop("There is no t-SNE transformation output files. Please locate it or generate a new one or use platform = 'R'")
                }else{
                  deleteTemp = FALSE
                }
              }
          if(deleteTemp){
            file.remove(paste(libpath,'/data.csv',sep=''))
            file.remove(paste(libpath,'/para.csv',sep=''))
            file.remove(paste(libpath,'/tsnetoutput.csv',sep=''))
            }
          }else{
            stop('t-SNE transform should be run in R or MATLAB')
          }
      }
    }
  cna<-vector()
  for(i in 1:ncol(tsner)){
    cna<-c(cna,paste('tsnet_',i,sep=''))
    }
  colnames(tsner)<-cna
  rownames(tsner)<-cells
  ttf.obj <- new(Class = "DimReduc",cell.embeddings = as.matrix(tsner))
  object@reductions[[reduction]]<-ttf.obj
  DefaultAssay(object@reductions[[reduction]])<-"RNA"
  return(object)
}

