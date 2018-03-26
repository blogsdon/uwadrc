loadData = function(loc = c('Updated TPP SNP Data.xlsx')){
  dataManifest <- list()
  dataManifest$genotypeAnnotations <- readxl::read_excel(loc[1],sheet = 'SNPs with ORs & Pathways')
  dataManifest$genotypes <- readxl::read_excel(loc[1],sheet = 'All Subjects')
  dataManifest$phenotypes <- readxl::read_excel(loc[1],sheet = 'Analysis')

  ####clean up snp annotation data
  dataManifest$genotypeAnnotations <- dplyr::select(dataManifest$genotypeAnnotations,
                                                    Chromosome,
                                                    Position,
                                                    `rs ID`,
                                                    `Closest Gene`,
                                                    Reference,
                                                    Alternate,
                                                    `Risk Allele`,
                                                    `Odds Ratio`,
                                                    `Risk Allele Frequency (In controls)`,
                                                    `Pathway 1`,
                                                    `Pathway 2`,
                                                    `Pathway 3`)
  colnames(dataManifest$genotypeAnnotations)[c(3,4,7:12)] <- c('rsID',
                                                               'Gene',
                                                               'RiskAllele',
                                                               'OddsRatio',
                                                               'RiskAlleleAlleleFrequency',
                                                               'PathwayAnnotation1',
                                                               'PathwayAnnotation2',
                                                               'PathwayAnnotation3')

  dataManifest$genotypeAnnotations$RiskAllele <- gsub('\\*','',dataManifest$genotypeAnnotations$RiskAllele)
  dataManifest$genotypes <- dplyr::select(dataManifest$genotypes,
                                          `Study ID`,
                                          `rs ID`,
                                          `Subject Genotype`,
                                          `Number of Risk Alleles`,
                                          `Odds Ratio`,
                                          `MAF`)
  colnames(dataManifest$genotypes) <- c('StudyID',
                                        'rsID',
                                        'Genotype',
                                        'nRiskAlleles',
                                        'OddsRatio',
                                        'MAF')
  dataManifest$genotypes <- dplyr::mutate(dataManifest$genotypes,
                                          logOddsRatio = log10(OddsRatio),
                                          sqrtMAF = sqrt(2*MAF*(1-MAF)))
  dataManifest$genotypes <- dplyr::mutate(dataManifest$genotypes,
                                          weight = logOddsRatio*sqrtMAF,
                                          genoWeight = weight*nRiskAlleles)
  dataManifest$genotypeSummary <- dplyr::group_by(dataManifest$genotypes,StudyID)
  dataManifest$genotypeSummary <- dplyr::summarise(dataManifest$genotypeSummary,
                                                   weightSum = sum(weight,na.rm=T),
                                                   genoWeightSum = sum(genoWeight,na.rm=T),
                                                   riskAlleleBurden = sum(nRiskAlleles,na.rm=T))
  dataManifest$genotypeSummary <- dplyr::mutate(dataManifest$genotypeSummary,
                                                GRS = 30*(genoWeightSum/weightSum))

  dataManifest$phenotypes <- dplyr::select(dataManifest$phenotypes,
                                           `Study ID`,
                                           `AOO`,
                                           `Diagnosis`,
                                           `Phenotype Notes`,
                                           `Research Variants`)
  colnames(dataManifest$phenotypes) <- c('StudyID',
                                         'AOO',
                                         'Diagnosis',
                                         'PhenotypeNotes',
                                         'ResearchVariants')
  dataManifest$individualSummary <- dplyr::left_join(dataManifest$genotypeSummary,
                                                     dataManifest$phenotypes,
                                                     by = 'StudyID')
  return(dataManifest)

}

