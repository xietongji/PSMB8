library(reshape2)
#step 1)Estimating TIME & 2)Transforming into binary variates
load("C:/Users/Lenovo/Desktop/Github/Biomarker selection.RData")
#step 3)~4) Wilcoxon rank sum test
CIBERSORTresults=CIBERSORTresults[CIBERSORTresults[,"P-value"]<0.05,]
CIBERSORTresults=CIBERSORTresults[,1:(ncol(CIBERSORTresults)-3)]
diff.results=data.frame()
for (Omics in c("SNV","CNAmp","CNDel","5mCMet","mRNA")) {
  ifelse(Omics=="SNV",(Groupingdata=SNVGrouping),
         ifelse(Omics=="CNAmp",(Groupingdata=CNVGrouping),
                ifelse(Omics=="CNDel",(Groupingdata=CNVGrouping),
                       ifelse(Omics=="5mCMet",(Groupingdata=MetGrouping),
                              ifelse(Omics=="mRNA",(Groupingdata=mRNAGrouping),NA)))))
  for (Algorithms in c("CIBERSORT","quanTIseq")) {
    ifelse(Algorithms=="CIBERSORT",(immunedata=CIBERSORTresults),
           ifelse(Algorithms=="quanTIseq",(immunedata=quanTIseqresults),NA))
    for(gene in c("PSMB8","PSMB9","ZBTB18")){
      print(paste0(Omics,Algorithms,gene))
      if(Omics=="CNAmp"){
        controlName=colnames(Groupingdata)[Groupingdata[gene,] %in% c(-2,-1,0)]
        researchName=colnames(Groupingdata)[Groupingdata[gene,] %in% c(1,2)]
      }else if(Omics=="CNDel"){
        controlName=colnames(Groupingdata)[Groupingdata[gene,] %in% c(2,1,0)]
        researchName=colnames(Groupingdata)[Groupingdata[gene,] %in% c(-1,-2)]
      }else{
        controlName=colnames(Groupingdata)[Groupingdata[gene,]==0]
        researchName=colnames(Groupingdata)[Groupingdata[gene,]==1]
      }
      controlName=intersect(rownames(immunedata),controlName)
      researchName=intersect(rownames(immunedata),researchName)
      subdiff.results=data.frame()
      if(as.numeric(length(researchName))==0){next}
      for(cell in colnames(immunedata)){
        wilcoxTest=wilcox.test(immunedata[controlName,cell],immunedata[researchName,cell])
        researchmedian=median(immunedata[researchName,cell])
        controlmedian=median(immunedata[controlName,cell])
        diff=researchmedian-controlmedian
        p=wilcoxTest$p.value
        subdiff.results=rbind(subdiff.results,
                              cbind(Omics=Omics,
                                    Algorithms=Algorithms,
                                    gene=gene,
                                    cell=cell,
                                    researchmedian=researchmedian,
                                    controlmedian=controlmedian,
                                    diff=diff,
                                    p=p))}
      pvalue=subdiff.results[,"p"]
      fdr=p.adjust(as.numeric(as.vector(pvalue)),method="fdr")
      subdiff.results=cbind(subdiff.results,fdr=fdr)
      diff.results=rbind(diff.results,subdiff.results)
    }
  }
}
#step 5)Statistic filtration & trend judgement
fourcellsResults=diff.results[diff.results$cell %in% c("T cells CD8",
                                                       "T.cells.CD8",
                                                       "T cells regulatory (Tregs)",
                                                       "Tregs",
                                                       "Macrophages M1",
                                                       "Macrophages.M1",
                                                       "Macrophages M2",
                                                       "Macrophages.M2"),]
fourcellsResults$cell=ifelse(fourcellsResults$cell %in% c("T cells CD8","T.cells.CD8"),"CD8T",
                             ifelse(fourcellsResults$cell %in% c("T cells regulatory (Tregs)","Tregs"),"Treg",
                                    ifelse(fourcellsResults$cell %in% c("Macrophages M1","Macrophages.M1"),"M1",
                                           ifelse(fourcellsResults$cell %in% c("Macrophages M2","Macrophages.M2"),"M2",NA))))
StatSigResults=fourcellsResults[(as.numeric(as.vector(fourcellsResults$p))<0.05 & as.numeric(as.vector(fourcellsResults$fdr))<0.25),]
LIST=list()
header=unique(paste0(fourcellsResults$Omics,"_",
                     fourcellsResults$cell,"_",
                     fourcellsResults$Algorithms))
for (setname in header) {
  Omics=strsplit(setname,"_")[[1]][1]
  cell=strsplit(setname,"_")[[1]][2]
  Algorithms=strsplit(setname,"_")[[1]][3]
  if(cell %in% c("CD8T","M1")){
    LIST[[paste0("TypeA_",setname)]]=StatSigResults[(StatSigResults$Omics==Omics &
                                                       StatSigResults$Algorithms==Algorithms &
                                                       StatSigResults$cell==cell &
                                                       as.numeric(as.vector(StatSigResults$diff))>0),"gene"]
    LIST[[paste0("TypeB_",setname)]]=StatSigResults[(StatSigResults$Omics==Omics &
                                                       StatSigResults$Algorithms==Algorithms &
                                                       StatSigResults$cell==cell &
                                                       as.numeric(as.vector(StatSigResults$diff))<0),"gene"]
  }else{
    LIST[[paste0("TypeC_",setname)]]=StatSigResults[(StatSigResults$Omics==Omics &
                                                       StatSigResults$Algorithms==Algorithms &
                                                       StatSigResults$cell==cell &
                                                       as.numeric(as.vector(StatSigResults$diff))>0),"gene"]
    LIST[[paste0("TypeD_",setname)]]=StatSigResults[(StatSigResults$Omics==Omics &
                                                       StatSigResults$Algorithms==Algorithms &
                                                       StatSigResults$cell==cell &
                                                       as.numeric(as.vector(StatSigResults$diff))<0),"gene"]
  }
}
#step 6) Consistent filtration of algorithms (between CIBERSORT & quanTIseq)  & trend (between CD8T & M1; Treg & M2)
LIST2=list()
header2=unique(paste0(fourcellsResults$Omics,"_",
                      fourcellsResults$cell))
for (setname in header2) {
  LIST2[[paste0("TypeA_",setname)]]=intersect(LIST[[paste0("TypeA_",setname,"_CIBERSORT")]],
                                              LIST[[paste0("TypeA_",setname,"_quanTIseq")]])
  LIST2[[paste0("TypeB_",setname)]]=intersect(LIST[[paste0("TypeB_",setname,"_CIBERSORT")]],
                                              LIST[[paste0("TypeB_",setname,"_quanTIseq")]])
  LIST2[[paste0("TypeC_",setname)]]=intersect(LIST[[paste0("TypeC_",setname,"_CIBERSORT")]],
                                              LIST[[paste0("TypeC_",setname,"_quanTIseq")]])
  LIST2[[paste0("TypeD_",setname)]]=intersect(LIST[[paste0("TypeD_",setname,"_CIBERSORT")]],
                                              LIST[[paste0("TypeD_",setname,"_quanTIseq")]])
}
LIST3=list()
header3=unique(fourcellsResults$Omics)
for (setname in header3) {
  LIST3[[paste0("TypeA_",setname)]]=intersect(LIST2[[paste0("TypeA_",setname,"_CD8T")]],
                                              LIST2[[paste0("TypeA_",setname,"_M1")]])
  LIST3[[paste0("TypeB_",setname)]]=intersect(LIST2[[paste0("TypeB_",setname,"_CD8T")]],
                                              LIST2[[paste0("TypeB_",setname,"_M1")]])
  LIST3[[paste0("TypeC_",setname)]]=intersect(LIST2[[paste0("TypeC_",setname,"_Treg")]],
                                              LIST2[[paste0("TypeC_",setname,"_M2")]])
  LIST3[[paste0("TypeD_",setname)]]=intersect(LIST2[[paste0("TypeD_",setname,"_Treg")]],
                                              LIST2[[paste0("TypeD_",setname,"_M2")]])
}
#step 7)~8) Consistent filtration of trend (between pro-tumor & anti-tumor immune cells)
fourcellsResultsnew=fourcellsResults
fourcellsResultsnew$diff=ifelse(as.numeric(as.vector(fourcellsResultsnew$p))<0.05 &
                                  as.numeric(as.vector(fourcellsResultsnew$fdr))<0.25,as.numeric(as.vector(fourcellsResultsnew$diff)),0)
fourcellsResultsCIBERSORTnew=dcast(fourcellsResultsnew[fourcellsResultsnew$Algorithms=="CIBERSORT",c("Omics","gene","cell","diff")],Omics+gene~cell)
fourcellsResultsquanTIseqnew=dcast(fourcellsResultsnew[fourcellsResultsnew$Algorithms=="quanTIseq",c("Omics","gene","cell","diff")],Omics+gene~cell)
fourcellsResultsCIBERSORTnew$Tcell=as.numeric(as.vector(fourcellsResultsCIBERSORTnew$CD8T))-as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Treg))
fourcellsResultsCIBERSORTnew$Macrophage=as.numeric(as.vector(fourcellsResultsCIBERSORTnew$M1))-as.numeric(as.vector(fourcellsResultsCIBERSORTnew$M2))
fourcellsResultsquanTIseqnew$Tcell=as.numeric(as.vector(fourcellsResultsquanTIseqnew$CD8T))-as.numeric(as.vector(fourcellsResultsquanTIseqnew$Treg))
fourcellsResultsquanTIseqnew$Macrophage=as.numeric(as.vector(fourcellsResultsquanTIseqnew$M1))-as.numeric(as.vector(fourcellsResultsquanTIseqnew$M2))
CIBERSORTLIST=list()
quanTIseqLIST=list()
for (setname in header3) {
  CIBERSORTLIST[[paste0("T+M+_",setname)]]=fourcellsResultsCIBERSORTnew[fourcellsResultsCIBERSORTnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsCIBERSORTnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Tcell))>0 &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Macrophage))>0,"gene"]
  CIBERSORTLIST[[paste0("T+M-_",setname)]]=fourcellsResultsCIBERSORTnew[fourcellsResultsCIBERSORTnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsCIBERSORTnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Tcell))>0 &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Macrophage))<0,"gene"]
  CIBERSORTLIST[[paste0("T-M+_",setname)]]=fourcellsResultsCIBERSORTnew[fourcellsResultsCIBERSORTnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsCIBERSORTnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Tcell))<0 &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Macrophage))>0,"gene"]
  CIBERSORTLIST[[paste0("T-M-_",setname)]]=fourcellsResultsCIBERSORTnew[fourcellsResultsCIBERSORTnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsCIBERSORTnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Tcell))<0 &
                                                                          as.numeric(as.vector(fourcellsResultsCIBERSORTnew$Macrophage))<0,"gene"]
  quanTIseqLIST[[paste0("T+M+_",setname)]]=fourcellsResultsquanTIseqnew[fourcellsResultsquanTIseqnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsquanTIseqnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Tcell))>0 &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Macrophage))>0,"gene"]
  quanTIseqLIST[[paste0("T+M-_",setname)]]=fourcellsResultsquanTIseqnew[fourcellsResultsquanTIseqnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsquanTIseqnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Tcell))>0 &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Macrophage))<0,"gene"]
  quanTIseqLIST[[paste0("T-M+_",setname)]]=fourcellsResultsquanTIseqnew[fourcellsResultsquanTIseqnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsquanTIseqnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Tcell))<0 &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Macrophage))>0,"gene"]
  quanTIseqLIST[[paste0("T-M-_",setname)]]=fourcellsResultsquanTIseqnew[fourcellsResultsquanTIseqnew$gene %in% c(LIST3[[paste0("TypeA_",setname)]],LIST3[[paste0("TypeB_",setname)]],LIST3[[paste0("TypeC_",setname)]],LIST3[[paste0("TypeD_",setname)]]) & 
                                                                          fourcellsResultsquanTIseqnew$Omics==setname &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Tcell))<0 &
                                                                          as.numeric(as.vector(fourcellsResultsquanTIseqnew$Macrophage))<0,"gene"]
}
LIST4=list()
for (setname in header3) {
  LIST4[[paste0("T+M+_",setname)]]=intersect(CIBERSORTLIST[[paste0("T+M+_",setname)]],
                                             quanTIseqLIST[[paste0("T+M+_",setname)]])
  LIST4[[paste0("T+M-_",setname)]]=intersect(CIBERSORTLIST[[paste0("T+M-_",setname)]],
                                             quanTIseqLIST[[paste0("T+M-_",setname)]])
  LIST4[[paste0("T-M+_",setname)]]=intersect(CIBERSORTLIST[[paste0("T-M+_",setname)]],
                                             quanTIseqLIST[[paste0("T-M+_",setname)]])
  LIST4[[paste0("T-M-_",setname)]]=intersect(CIBERSORTLIST[[paste0("T-M-_",setname)]],
                                             quanTIseqLIST[[paste0("T-M-_",setname)]])
}
#step 9) Multi-omics biomarkers
intersect(LIST4[["T+M+_mRNA"]],LIST4[["T-M-_5mCMet"]])
intersect(LIST4[["T+M+_5mCMet"]],LIST4[["T-M-_mRNA"]])
