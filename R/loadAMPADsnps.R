loadAMPADsnps <- function(){
  foo <- data.table::fread('ROSMAP_Mayo_Extracted_SNPs_forBen.txt',
                           data.table=F)
  bar <- data.table::fread('snps',
                           data.table=F)
  colnames(bar) <- c("rsId",
                     "loc",
                     "study",
                     "Ref",
                     "Alt",
                     "Risk",
                     "OR",
                     "Gene",
                     "Pathway")

  toFlip <- bar$rsId[which(bar$Ref == bar$Risk)]

  fooNew <- foo
  fooNew[,toFlip] <- 2-foo[,toFlip]
  #convert imputed genotypes to scores

  ors <- bar$OR
  names(ors) <- bar$rsId
  scores <- as.matrix(fooNew[,names(ors)])%*%ors
  ind <- which(bar$Pathway=='endosomal')
  endoScores <- as.matrix(fooNew[,names(ors)[ind]])%*%ors[ind]

  ind <- which(bar$Pathway=='immune')
  immuneScores <- as.matrix(fooNew[,names(ors)[ind]])%*%ors[ind]

  ind <- which(bar$Pathway=='neuronal')
  neuroScores <- as.matrix(fooNew[,names(ors)[ind]])%*%ors[ind]

  scoreDf <- data.frame(riskScore = scores,
                        endoScore = endoScores,
                        immuneScore = immuneScores,
                        neuroScore = neuroScores,
                        ID = fooNew$ID)

  scoreDf<-dplyr::left_join(fooNew,scoreDf)

  return(scoreDf)
}
