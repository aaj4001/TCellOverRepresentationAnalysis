library(openxlsx)
library(clusterProfiler)
library(readr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(viridis)

########################################################################################################################
## Imports Gene Sets (with Universe)

Arrange_ForScoring <- function(SignatureList,AllGenes){
  SigList_ForScoring = list()
  for(i in 1:length(SignatureList)) SigList_ForScoring[[i]] = list(Signature = SignatureList[[i]],AllGenes = AllGenes)
  names(SigList_ForScoring) = names(SignatureList)
  
  SigList_ForScoring
}

Hacohen_DEGs = msigdbi::read.gmt("Hacohen_StateDEGs_Human.gmt")$genesets
HacohenAllGenes = clusterProfiler::read.gmt("Hacohen_AllGenes.gmt")

Hacohen_DEGs = Arrange_ForScoring(Hacohen_DEGs,HacohenAllGenes)

## Imports Query Sets Figure 4A
MurineSigs = msigdbi::read.gmt("MurineSigs.gmt")$genesets
names(MurineSigs)[c(7,8)] = c("Act_Dys","Naive_Mem")

YellowFever = msigdbi::read.gmt("HumanViralTCellsSigs.gmt")$genesets[c(3,4,7,8)]
               

########################################################################################################################
## Overrepresentation Functions 

ListtoPathway = function(InputList,AllGenes){
  term = NULL
  gene = NULL
  for(i in 1:length(InputList)) {
    term = c(term,rep(names(InputList)[i],length(InputList[[i]])))
    gene = c(gene,InputList[[i]])
  }
  
  PathwayOutput = data.frame(term,gene);rownames(PathwayOutput) = NULL
  PathwayOutput$gene = as.character(PathwayOutput$gene)
  
  PathwayOutput = rbind(PathwayOutput,AllGenes)
  
  PathwayOutput
}

#InputSigList = Hacohen_DEGs;InputQueryList = YellowFever;i = 1
## Overrepresentation Calculation
OverRepresent_Results <- function(InputSigList, InputQueryList){
  OverRepresent_Result = list()
  for(i in 1:length(InputSigList)){
    OverRepresent_Result[[i]] = enricher(InputSigList[[i]]$Signature,
                                         TERM2GENE = ListtoPathway(InputQueryList,InputSigList[[i]]$AllGenes),
                                         minGSSize = 1,maxGSSize = 10000,pvalueCutoff = 2,qvalueCutoff = 2)
  }
  
  names(OverRepresent_Result) = names(InputSigList)
  OverRepresent_Result
}

OverRepresent_Overlap <- function(InputSigList,InputQueryList){
  MergedResults = matrix(nrow = length(InputSigList),ncol = length(InputQueryList))
  OverRepresent_Result = OverRepresent_Results(InputSigList,InputQueryList)
  
  for(i in 1:length(InputSigList)){
    if(!is.null(OverRepresent_Result[[i]])){
      GeneRatio = OverRepresent_Result[[i]]@result[,c(1,9)]
      Result = merge(data.frame(PathwayName = names(InputQueryList),sort = 1:length(InputQueryList)),
                     GeneRatio,by = 1,all.x = TRUE)
      Result = Result[order(Result$sort),]
      MergedResults[i,] = as.character(Result$Count)
    }
  }
  
  rownames(MergedResults) = names(InputSigList);colnames(MergedResults) = names(InputQueryList)
  MergedResults[is.na(MergedResults)] = ""
  
  MergedResults
}

OverRepresent_OverlapGenes <- function(InputSigList,InputQueryList){
  MergedResults = matrix(nrow = length(InputSigList),ncol = length(InputQueryList))
  OverRepresent_Result = OverRepresent_Results(InputSigList,InputQueryList)
  
  for(i in 1:length(InputSigList)){
    if(!is.null(OverRepresent_Result[[i]])){
      GeneOverlaps = OverRepresent_Result[[i]]@result[,c(1,8)]
      Result = merge(data.frame(PathwayName = names(InputQueryList),sort = 1:length(InputQueryList)),
                     GeneOverlaps,by = 1,all.x = TRUE)
      Result = Result[order(Result$sort),]
      MergedResults[i,] = Result$geneID
    }
  }
  
  rownames(MergedResults) = names(InputSigList);colnames(MergedResults) = names(InputQueryList)
  MergedResults[is.na(MergedResults)] = "N/A"
  
  MergedResults
}

OverRepresent_PVal <- function(InputSigList,InputQueryList){
  MergedResults = matrix(nrow = length(InputSigList),ncol = length(InputQueryList))
  OverRepresent_Result = OverRepresent_Results(InputSigList,InputQueryList)
  
  for(i in 1:length(InputSigList)){
    if(!is.null(OverRepresent_Result[[i]])){
      Pvalue = OverRepresent_Result[[i]]@result[,c(1,5)]
      Result = merge(data.frame(PathwayName = names(InputQueryList),sort = 1:length(InputQueryList)),
                     Pvalue,by = 1,all.x = TRUE)
      Result = Result[order(Result$sort),]
      MergedResults[i,] = Result$pvalue
    }
  }
  
  rownames(MergedResults) = names(InputSigList);colnames(MergedResults) = names(InputQueryList)
  MergedResults[is.na(MergedResults)] = 1.0
  
  MergedResults
}

OverRepresent_OddsRatio <- function(InputSigList,InputQueryList){
  MergedResults = matrix(nrow = length(InputSigList),ncol = length(InputQueryList))
  OverRepresent_Result = OverRepresent_Results(InputSigList,InputQueryList)
  
  for(i in 1:length(InputSigList)){
    if(!is.null(OverRepresent_Result[[i]])){
      GeneRatio = OverRepresent_Result[[i]]@result[,c(1,3)]
      GeneRatio$GeneRatio = sapply(strsplit(GeneRatio$GeneRatio,split = "/",fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2]))
      
      BgRatio = OverRepresent_Result[[i]]@result[,c(1,4)]
      BgRatio$BgRatio = sapply(strsplit(BgRatio$BgRatio,split = "/",fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2]))
      
      OddsRatio = merge(GeneRatio,BgRatio,by = 1)
      OddsRatio$OddsRatio = OddsRatio$GeneRatio/OddsRatio$BgRatio
      
      Result = merge(data.frame(PathwayName = names(InputQueryList),sort = 1:length(InputQueryList)),
                     OddsRatio,by = 1,all.x = TRUE)
      Result = Result[order(Result$sort),]
      MergedResults[i,] = Result$OddsRatio
    }
  }  
  
  rownames(MergedResults) = names(InputSigList);colnames(MergedResults) = names(InputQueryList)
  MergedResults[is.na(MergedResults)] = 0
  MergedResults
}

# InputSigList = TOX_Sig[rev(1:2)];InputQueryList = MurineSigs;UpperOddsRatio = 6.1
PicnicPlot <- function(InputSigList,InputQueryList,FDRCutoff = 0.05,
                       UpperPVal = NULL, UpperOddsRatio = NULL,
                       OverlapNumbers = TRUE,FlipVertical = FALSE){
  InputFDRs = data.frame(OverRepresent_PVal(InputSigList,InputQueryList))
  InputOddsRatios = data.frame(OverRepresent_OddsRatio(InputSigList,InputQueryList))
  InputOddsRatios = log2(InputOddsRatios)
  InputResults = InputFDRs
  InputResults = -1*log10(InputResults)
  InputOverlaps = data.frame(OverRepresent_Overlap(InputSigList,InputQueryList),stringsAsFactors = FALSE)
  
  ## FDR adjusts p value (groups by row)
  for(i in 1:nrow(InputFDRs)) InputFDRs[i,] = p.adjust(InputFDRs[i,],method = "BH")
  
  Ordering = names(InputSigList)
  
  ## Merges both results together
  InputResults$TCellState_Sigs = rownames(InputFDRs)
  InputFDRs$TCellState_Sigs = rownames(InputFDRs)
  InputOddsRatios$TCellState_Sigs = rownames(InputFDRs)
  InputOverlaps$TCellState_Sigs = rownames(InputFDRs)
  
  
  InputResults = reshape2::melt(InputResults,id.vars = "TCellState_Sigs",variable.name = "Test_Signatures",value.name = "neglog10_pval")
  InputFDRs = reshape2::melt(InputFDRs,id.vars = "TCellState_Sigs",variable.name = "Test_Signatures",value.name = "p.adj")
  InputOddsRatios = reshape2::melt(InputOddsRatios,id.vars = "TCellState_Sigs",variable.name = "Test_Signatures",value.name = "log2_OddsRatio")
  InputOverlaps = reshape2::melt(InputOverlaps,id.vars = "TCellState_Sigs",variable.name = "Test_Signatures",value.name = "Overlap")
  
  InputResults$Merger = paste(InputResults$TCellState_Sigs,InputResults$Test_Signatures,sep = "_")
  InputFDRs$Merger = paste(InputFDRs$TCellState_Sigs,InputFDRs$Test_Signatures,sep = "_");InputFDRs = InputFDRs[,-(1:2)]
  InputOddsRatios$Merger = paste(InputOddsRatios$TCellState_Sigs,InputOddsRatios$Test_Signatures,sep = "_");InputOddsRatios = InputOddsRatios[,-(1:2)]
  InputOverlaps$Merger = paste(InputOverlaps$TCellState_Sigs,InputOverlaps$Test_Signatures,sep = "_");InputOverlaps = InputOverlaps[,-(1:2)]
  
  InputResults = merge(InputResults,InputFDRs,by.x = 4,by.y = 2)
  InputResults = merge(InputResults,InputOddsRatios,by.x = 1,by.y = 2)
  InputResults = merge(InputResults,InputOverlaps,by.x = 1,by.y = 2)[,-1]
  
  InputResults$TCellState_Sigs = factor(InputResults$TCellState_Sigs,levels = Ordering)
  InputResults$Test_Signatures = factor(InputResults$Test_Signatures, levels= names(InputQueryList))
  InputResults$Sig = "ns"
  InputResults$Sig[InputResults$p.adj < FDRCutoff] = "sig"

  PicPlot = ggplot(InputResults,aes(x = Test_Signatures,y = TCellState_Sigs,color = log2_OddsRatio,fill = log2_OddsRatio,size = neglog10_pval))

  if(is.null(UpperOddsRatio)&is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_color_viridis(option = "plasma",direction = -1) +
      scale_fill_viridis(option = "plasma",direction = -1) +
      geom_point(size = 0)
  } else if (!is.null(UpperOddsRatio)&is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_color_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      scale_fill_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      geom_point(size = 0)
  } else if(is.null(UpperOddsRatio)&!is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,UpperPVal)) +
      scale_color_viridis(option = "plasma",direction = -1) +
      scale_fill_viridis(option = "plasma",direction = -1) +
      geom_point(size = 0)
  } else if(!is.null(UpperOddsRatio)&!is.null(UpperPVal)){
    PicPlot = PicPlot +
      scale_size_continuous(limits = c(0,UpperPVal)) +
      scale_color_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      scale_fill_viridis(option = "plasma",direction = -1,limits = c(0,UpperOddsRatio)) +
      geom_point(size = 0)
  }  
  LegPlot = PicPlot
  PicPlot = PicPlot + 
    geom_point(data = InputResults[InputResults$Sig == "ns",],color = "grey",fill = "grey",shape = 22) +
    geom_point(data = InputResults[InputResults$Sig == "sig",],shape = 21) +
    theme_bw() + 
    # scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90,vjust =0.5,hjust = 1),axis.title = element_blank())
    
  LegPlot = LegPlot + 
    geom_point() + 
    theme_bw() + 
    # scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),axis.title = element_blank())
  
  if(OverlapNumbers){
    PicPlot = PicPlot + geom_text(aes(label = Overlap),data = InputResults[InputResults$Sig=="sig",],size = 2,color = "white")
    LegPlot = LegPlot + geom_text(aes(label = Overlap),data = InputResults[InputResults$Sig=="sig",],size = 2,color = "white")
  }
  
  if(FlipVertical){
    PicPlot = PicPlot + coord_flip() + 
      scale_y_discrete(limits = rev(levels(InputResults$TCellState_Sigs)),position = "right") +
      theme(axis.text.x = element_text(hjust = 0,vjust = 0.5)) +
      scale_x_discrete(limits = rev(levels(InputResults$Test_Signatures)),position = "top")
      
    LegPlot = LegPlot + coord_flip() + 
      scale_y_discrete(limits = rev(levels(InputResults$TCellState_Sigs)),position = "right") +
      theme(axis.text.x = element_text(hjust = 0,vjust = 0.5)) +
      scale_x_discrete(limits = rev(levels(InputResults$Test_Signatures)),position = "top")
  }
  
  Plots = list(PicnicPlot = PicPlot,LegendPlot = LegPlot)
  Plots
}


########################################################################################################################
## Overrepresentation Analysis Results

max(log2(OverRepresent_OddsRatio(Hacohen_DEGs,YellowFever)))
max(-log10(OverRepresent_PVal(Hacohen_DEGs,YellowFever)))

HacohenDEG_PicPlot_Murine = PicnicPlot(Hacohen_DEGs[6:1],MurineSigs,FlipVertical = TRUE,
                                       UpperPVal = 36,UpperOddsRatio = 5.2)


HacohenDEG_PicPlot_YF = PicnicPlot(Hacohen_DEGs[6:1],YellowFever,FlipVertical = TRUE,
                                   UpperOddsRatio = 4.7,UpperPVal = 61)
HacohenDEG_PicPlot_YF$PicnicPlot

cowplot::plot_grid(HacohenDEG_PicPlot_Murine$PicnicPlot,
                   HacohenDEG_PicPlot_YF$PicnicPlot + theme(axis.text.x = element_blank()),
                   nrow = 2,align = "v",rel_heights = c(1,0.2))



Scale = 0.7

pdf("4AB_OverrepresentionAnalysis_Hacohen_v4.pdf",width = 6.9*Scale,height = 10.45*Scale)
cowplot::plot_grid(HacohenDEG_PicPlot_Murine$PicnicPlot,
                   HacohenDEG_PicPlot_YF$PicnicPlot + theme(axis.text.x = element_blank()),
                   nrow = 2,align = "v",rel_heights = c(1,0.2))
dev.off()

pdf("4AB_OverrepresentionAnalysis_Hacohen_v4Legend.pdf",width = 6.9*Scale,height = 30*Scale)
cowplot::plot_grid(HacohenDEG_PicPlot_Murine$LegendPlot,
                   HacohenDEG_PicPlot_YF$LegendPlot + theme(axis.text.x = element_blank()),
                   nrow = 2,align = "v",rel_heights = c(1,0.2))
dev.off()
