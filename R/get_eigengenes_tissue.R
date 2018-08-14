get_eigengenes_tissue <- function(synId = 'syn11932957',tissue = 'DLPFC'){
  #aggMods <- rSynapseUtilities::loadFullTable(synId)
  #synapseClient::synapseLogin()
  synapser::synLogin()
  #mods <- synapseClient::synTableQuery(paste0("select * from ",synId," where brainRegion = \'",tissue,'\''))@values
  mods <- synapser::synTableQuery(paste0("select * from ",synId," where brainRegion = \'",tissue,'\''))$asDataFrame()
  mods <- mods[,-c(1,2)]
  #library(AMPAD)
  #geneExpressionForAnalysis <- AMPAD::pullExpressionAndPhenoWinsorized()
  load('gexpr.rda')
  names(geneExpressionForAnalysis) <- c('TCX','CBE','DLPFC','FP','STG','PHG','IFG')
  computeEigengene <- function(br,geneExp,moduleDefinitions){

    geneExp <- geneExp[[br]]

    #get modules
    #mods <- dplyr::filter(moduleDefinitions,brainRegion==br)
    #convert modules into list of genes
    modsDefs <- lapply(unique(mods$ModuleNameFull),
                       utilityFunctions::listify,
                       mods$GeneID,
                       mods$ModuleNameFull)

    names(modsDefs) <- unique(mods$ModuleNameFull)

    internal <- function(mod,modsDefs,geneExp){
      geneExpMod <- dplyr::select(geneExp,modsDefs[[mod]])
      geneExpMod <- scale(geneExpMod)
      foo <- svd(geneExpMod)
      eigenGenes <- as.matrix(foo$u[,1:15])
      eigenValues <- foo$d
      weights <- foo$v
      rownames(weights) <- colnames(geneExpMod)
      colnames(eigenGenes) <- paste0('pc',1:15)
      rownames(eigenGenes) <- geneExp$aSampleId
      #res <- cor(eigenGenes,geneExpMod)
      return(list(eigenGenes=eigenGenes,
                  eigenValues=eigenValues,
                  weights=weights))
    }

    full_res<-lapply(names(modsDefs),internal,modsDefs,geneExp)
    names(full_res) <- names(modsDefs)
    return(full_res)
  }

  fullList<-computeEigengene(tissue,geneExpressionForAnalysis,mods)
  #fullList<-lapply(names(geneExpressionForAnalysis)[keep],computeEigengene,geneExpressionForAnalysis,aggMods)
  #names(fullList) <- names(geneExpressionForAnalysis)
  return(fullList)
}
