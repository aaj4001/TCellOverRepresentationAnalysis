library(openxlsx)
library(clusterProfiler)
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

########################################################################################################################
## Imports Gene Sets (with Universe)

Hacohen_DEGs = msigdbi::read.gmt("Hacohen_StateDEGs_Human.gmt")$genesets
HacohenAllGenes = msigdbi::read.gmt("Hacohen_AllGenes.gmt")$genesets

## Imports Query Sets Figure 4A
MurineSigs = msigdbi::read.gmt("MurineSigs.gmt")$genesets
names(MurineSigs)[c(7,8)] = c("Act_Dys","Naive_Mem")

YellowFever = msigdbi::read.gmt("HumanViralTCellsSigs.gmt")$genesets[c(3,4,7,8)]
     
Test = read.gmt("MurineSigs.gmt")          
########################################################################################################################
## Overrepresentation Functions 


OverRepresent_Results <- function(QueryDEGs, SignatureList, Universe){
  ## Arranges SignatureList and Universe in right format for TERM2GENE
  SigList = lapply(seq_along(SignatureList),function(i, Signatures, Names) {
    data.frame(term = Names[i], gene = SignatureList[[i]])},
    Signatures = SignatureList, Names = names(SignatureList)) %>%
    bind_rows() %>% rbind(data.frame(term = "GeneUniverse",gene = Universe[[1]])) %>%
    mutate(term = factor(term, levels = c(names(SignatureList), "GeneUniverse")))
  
  ## Runs Overrepresentation Analysis
  OverRep_Result = lapply(QueryDEGs,enricher,TERM2GENE = SigList, minGSSize = 1,
                          maxGSSize = Inf,pvalueCutoff = 2,qvalueCutoff = 2)
  
  ## Extracts Results and adds Odds Ratio
  OverRep_Result = lapply(OverRep_Result, function(x) x@result)
  OverRep_Result = lapply(OverRep_Result, function(x) x %>% 
                            mutate(OddsRatio = sapply(strsplit(GeneRatio,split = "/",fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2]))/
                                     sapply(strsplit(BgRatio,split = "/",fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2]))))
  OverRep_Result
}



OverRepresent_Plot <- function(QueryDEGs, SignatureList, Universe,
                               FDRCutoff = 0.05, MakeSquare = F,
                               UpperPVal = NULL, UpperOddsRatio = NULL,
                               OverlapNumbers = TRUE,FlipVertical = FALSE){

  OverRep_Result = OverRepresent_Results(QueryDEGs, SignatureList, Universe)
  
  ## Adds in pathways with zero overlap, sets p value for non-overlapping to 1, Odds Ratio to 0 and calculates FDR
  OverRep_Result %<>% lapply(merge,y = names(SignatureList), by = 1, all.y = T)
  OverRep_Result %<>% lapply(function(x) {
    x$pvalue[is.na(x$Description)] = 1
    x$OddsRatio[is.na(x$Description)] = 0
    x %>% mutate(p.adjust = p.adjust(pvalue, method = "BH"),
                 ID = factor(ID, names(SignatureList)))
  })
  
  OverRep_Result = lapply(seq_along(OverRep_Result),
                          function(i, List, Names) List[[i]] %>% 
                            mutate(QueryDEG = Names[i]), 
                          List = OverRep_Result, Names = names(OverRep_Result)) %>%
    bind_rows() %>% mutate(QueryDEG = factor(QueryDEG, names(QueryDEGs)))
  
  ## Plots Enrichment (picnic plot)
  
  PicPlot = ggplot(OverRep_Result,aes(x = ID,y = QueryDEG,
                                      color = log2(OddsRatio),fill = log2(OddsRatio),size = -log10(pvalue)))

  if(is.null(UpperOddsRatio)&is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,NA)) +
      scale_color_viridis(option = "plasma",direction = -1, limits = c(0,NA)) +
      scale_fill_viridis(option = "plasma",direction = -1, limits = c(0,NA)) +
      geom_point(size = 0)
  } else if (!is.null(UpperOddsRatio)&is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,NA)) +
      scale_color_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      scale_fill_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      geom_point(size = 0)
  } else if(is.null(UpperOddsRatio)&!is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,UpperPVal)) +
      scale_color_viridis(option = "plasma",direction = -1, limits = c(0,NA)) +
      scale_fill_viridis(option = "plasma",direction = -1, limits = c(0,NA)) +
      geom_point(size = 0)
  } else if(!is.null(UpperOddsRatio)&!is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,UpperPVal)) +
      scale_color_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      scale_fill_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      geom_point(size = 0)
  }  
  
  PicPlot = PicPlot + 
    geom_point(data = OverRep_Result %>% subset(p.adjust > FDRCutoff),color = "grey",fill = "grey",shape = 22,show.legend = F) +
    geom_point(data = OverRep_Result %>% subset(p.adjust <= FDRCutoff)) +
    theme_bw() + theme(axis.title = element_blank()) 
  
  
  if(OverlapNumbers){
    PicPlot = PicPlot + geom_text(aes(label = Count),data = OverRep_Result %>% subset(p.adjust <= FDRCutoff),size = 2,color = "white")
  }
  
  if(FlipVertical){
    PicPlot = PicPlot + scale_x_discrete(limits = rev) +
      scale_y_discrete(position = "right") + 
      coord_flip() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
            legend.box = "horizontal")
  } else{
    PicPlot = PicPlot + scale_y_discrete(limits = rev) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.box = "vertical", legend.position = "top")
  }
  
  if(MakeSquare & !FlipVertical) {PicPlot = PicPlot + theme(aspect.ratio = length(QueryDEGs)/length(SignatureList))
  } else if(MakeSquare & FlipVertical) PicPlot = PicPlot + theme(aspect.ratio = length(SignatureList)/length(QueryDEGs))

  PicPlot
}

OverRepresent_MultiPlot <- function(QueryDEGs, Listof_SignatureLists, Universe,
                                    FDRCutoff = 0.05, MakeSquare = F, Silent = F,
                                    UpperPVal = NULL, UpperOddsRatio = NULL,
                                    OverlapNumbers = TRUE,FlipVertical = FALSE){
  OverRep_Results = lapply(Listof_SignatureLists, OverRepresent_Results, QueryDEGs = QueryDEGs, Universe = Universe)

  OverRep_Results %<>% unlist(recursive = F) %>% bind_rows()

  if(is.null(UpperPVal)) UpperPVal = max(-log10(OverRep_Results$pvalue))
  if(is.null(UpperOddsRatio)) UpperOddsRatio = max(log2(OverRep_Results$OddsRatio))
  Widths = sapply(Listof_SignatureLists, length)
  
  MultiPlot = lapply(Listof_SignatureLists, OverRepresent_Plot, QueryDEGs = QueryDEGs, Universe = Universe,
                     FDRCutoff = FDRCutoff, MakeSquare = MakeSquare, UpperPVal = UpperPVal, UpperOddsRatio = UpperOddsRatio,
                     OverlapNumbers = OverlapNumbers, FlipVertical = FlipVertical)
  
  if(!Silent){
    if(!FlipVertical) {
      for(i in 2:length(MultiPlot)) MultiPlot[[i]] = MultiPlot[[i]] + theme(axis.text.y = element_blank())
      plot_grid(plotlist = MultiPlot, rel_widths = Widths, align = "h",axis = "lb")
    } else {
      for(i in 2:length(MultiPlot)) MultiPlot[[i]] = MultiPlot[[i]] + theme(axis.text.x = element_blank())
      plot_grid(plotlist = MultiPlot, nrow = length(MultiPlot),rel_heights = Widths, align = "v",axis = "lb")
    }
  } else{
    MultiPlot
  }
}

########################################################################################################################
## Overrepresentation Analysis Use Case

OverRepresent_Plot(Hacohen_DEGs, YellowFever, HacohenAllGenes)

OverRepresent_MultiPlot(Hacohen_DEGs, list(MurineSigs, YellowFever), HacohenAllGenes)
