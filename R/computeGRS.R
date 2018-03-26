computeGRS <- function(x,y){
  gs <- dplyr::group_by(x,StudyID)
  gs <- dplyr::summarise(gs,
                                                   weightSum = sum(weight,na.rm=T),
                                                   genoWeightSum = sum(genoWeight,na.rm=T),
                                                   riskAlleleBurden = sum(nRiskAlleles,na.rm=T))
  gs <- dplyr::mutate(gs,
                                                GRS = 30*(genoWeightSum/weightSum))
  colnames(gs)[2:5] <- paste0(colnames(gs[2:5]),y)
  return(gs)
}
