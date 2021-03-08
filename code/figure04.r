
rm(list=ls())
library(egg)
library(gridExtra)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(pracma)
library(scales)
require(GGally)
require(gtools)

# set the working directory to the exact code base location ------------------
setwd(paste0(path.expand("~"),"/Desktop/shaikh-sears-clin-cancer-research-2021/code/"))
# ------------------
analysisdir = "../analysis/figure04/"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ------------------
pdfdir     = analysisdir
# ------------------
coefficients <- list(
  "Akkermansia_muciniphila"=1,
  "Alistipes_indistinctus"=1,
  "Anaerostipes_hadrus"=1,
  "Anaerotruncus_colihominis"=-1,
  "Bacteroides_caccae"=1,
  "Bacteroides_coprocola"=-1,
  "Bacteroides_fragilis"=-1,
  "Bacteroides_thetaiotaomicron"=-1,
  "Bacteroides_uniformis"=-0.1,
  "Bifidobacterium_adolescentis"=1,
  "Bifidobacterium_longum"=1,
  "Clostridium_hathewayi"=-1,
  "Clostridium_hylemonae"=-1,
  "Clostridium_methylpentosum"=-1,
  "Collinsella_aerofaciens"=1,
  "Dorea_formicigenerans"=1,
  "Enterococcus_faecium"=1,
  "Faecalibacterium_prausnitzii"=1,
  "g.Akkermansia"=1,
  "g.Bifidobacterium"=1,
  "g.Enterococcus"=1,
  "g.Eubacterium"=1,
  "g.Faecalibacterium"=1,
  "g.Lactobacillus"=1,
  "g.Parabacteroides"=-1,
  "g.Ruminococcus"=1,
  "Gemmiger_formicilis"=1,
  "Holdemania_filiformis"=1,
  "Klebsiella_pneumoniae"=1,
  "Megasphaera_micronuciformis"=-1,
  "Methanobrevibacter_smithii"=1,
  "Oribacterium_sinus"=-1,
  "Oxalobacter_formigenes"=1,
  "Parabacteroides_merdae"=-1,
  "Parasutterella_excrementihominis"=-1,
  "Prevotella_buccalis"=1,
  "Roseburia_hominis"=1,
  "Roseburia_intestinalis"=-1,
  "Ruminococcus_obeum"=-1,
  "Scardovia_wiggsiae"=-1,
  "Streptococcus_parasanguinis"=1,
  "Veillonella_parvula"=-1
)

indexes = list(
  "lit-r"   = c("g.Faecalibacterium","Akkermansia_muciniphila","Gemmiger_formicilis","Bacteroides_caccae","Streptococcus_parasanguinis","Holdemania_filiformis","Bacteroides_thetaiotaomicron","Bifidobacterium_longum","Collinsella_aerofaciens","Enterococcus_faecium"),
  "lit-nr"  = c("g.Parabacteroides", "Ruminococcus_obeum", "Roseburia_intestinalis"),
  "lit-agg" = c("g.Faecalibacterium","Akkermansia_muciniphila", "Gemmiger_formicilis","Bacteroides_caccae","Streptococcus_parasanguinis","Holdemania_filiformis","Bacteroides_thetaiotaomicron","Bifidobacterium_longum","Collinsella_aerofaciens","Enterococcus_faecium","g.Parabacteroides", "Ruminococcus_obeum", "Roseburia_intestinalis"),

  "lit-melanoma-r"    = c("g.Faecalibacterium","Gemmiger_formicilis", "Bacteroides_caccae","Streptococcus_parasanguinis","Holdemania_filiformis","Bacteroides_thetaiotaomicron","Bifidobacterium_longum","Collinsella_aerofaciens","Enterococcus_faecium"),
  "lit-melanoma-nr"   = c("Ruminococcus_obeum", "Roseburia_intestinalis"),
  "lit-melanoma-agg"  = c("g.Faecalibacterium","Gemmiger_formicilis", "Bacteroides_caccae","Streptococcus_parasanguinis","Holdemania_filiformis","Bacteroides_thetaiotaomicron","Bifidobacterium_longum","Collinsella_aerofaciens","Enterococcus_faecium", "Ruminococcus_obeum", "Roseburia_intestinalis"),

  "re-analysis-r"     = c(
    "Faecalibacterium_prausnitzii",
    "Akkermansia_muciniphila",
    "Prevotella_buccalis",
    "Roseburia_hominis",
    "Oxalobacter_formigenes"
  ),
  "re-analysis-nr"    =
  c(
    "Veillonella_parvula",
    "Clostridium_hathewayi",
    "Clostridium_hylemonae",
    "Bacteroides_coprocola",
    "Bacteroides_fragilis",
    "Parasutterella_excrementihominis",
    "Clostridium_methylpentosum",
    "Scardovia_wiggsiae",
    "Oribacterium_sinus",
    "Bacteroides_thetaiotaomicron",
    "Megasphaera_micronuciformis",
    "Bacteroides_uniformis"
  ),
  "re-analysis-agg"     = c(
    "Faecalibacterium_prausnitzii",
    "Akkermansia_muciniphila",
    "Prevotella_buccalis",
    "Roseburia_hominis",
    "Oxalobacter_formigenes",
    "Bacteroides_coprocola",
    "Bacteroides_fragilis",
    "Bacteroides_thetaiotaomicron",
    "Bacteroides_uniformis",
    "Clostridium_hathewayi",
    "Clostridium_hylemonae",
    "Clostridium_methylpentosum",
    "Megasphaera_micronuciformis",
    "Oribacterium_sinus",
    "Parasutterella_excrementihominis",
    "Scardovia_wiggsiae",
    "Veillonella_parvula"
  )
)

tablefile = "../analysis/A03-compile-indicator-taxa/out.A03.indicator-taxa.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, as.is=TRUE)
A <- A[A$UseforMeta=="Yes",]
A$Study = factor(A$Study)

indexData = c()
for (ci in names(indexes)){
  for (j in 1:nrow(A)){
    jrow       = c()
    indexvalue = 0

    if (!grepl("log.reg", ci)){
      for (ti in indexes[[ci]]){
        indexvalue = indexvalue + (coefficients[[ti]]*as.numeric(A[j,ti]))
      }
    }else{
      indexvalue = A[j,ci]
    }
    jrow = c(as.character(A[j,1]), as.character(A[j,2]), as.character(A[j,4]), as.character(A[j,7]), ci, indexvalue)
    indexData = rbind(indexData, jrow)
  }
}
colnames(indexData) = c("SampleID", "Study", "CancerType", "Response", "Type", "Value")
indexData           = data.frame(indexData)
indexData$Value     = as.numeric(as.character(indexData$Value))
write.table(indexData, file=paste(analysisdir, "index-results.txt", sep=""), col.names=TRUE, row.names=FALSE, sep="\t")

# -----------------------------------
iresultsT = c()
for (type in levels(indexData$Type))
for (study in levels(indexData$Study)){
  Set1 = indexData[indexData$Type==type & indexData$Study==study, ]
  lBs  = c(min(Set1$Value)-1, unique(sort(Set1$Value)), max(Set1$Value)+1)
  for (i in 1:length(lBs)){
    # now lBs[i] is our potential cutpoint to distinguish above sc or below sc
    Set1indexDataboveThresh <- Set1[Set1$Value >= lBs[i], ]$Response
    Set1BelowThresh         <- Set1[Set1$Value <  lBs[i], ]$Response

    # construct a 2x2 contigency table based on this information
    twobytwo      = rbind(c(0,0),c(0,0))
    twobytwo[1,1] = sum(Set1indexDataboveThresh == "R")
    twobytwo[1,2] = sum(Set1indexDataboveThresh == "NR")
    twobytwo[2,1] = sum(Set1BelowThresh == "R")
    twobytwo[2,2] = sum(Set1BelowThresh == "NR")

    SN  = 100*twobytwo[1,1]/(twobytwo[1,1]+twobytwo[2,1])
    SP  = 100*twobytwo[2,2]/(twobytwo[1,2]+twobytwo[2,2])
    PPV = 100*twobytwo[1,1]/(twobytwo[1,1]+twobytwo[1,2])
    NPV = 100*twobytwo[2,2]/(twobytwo[2,1]+twobytwo[2,2])
    ACC = 100*(twobytwo[1,1]+twobytwo[2,2]) / (twobytwo[1,1]+twobytwo[2,1]+twobytwo[1,2]+twobytwo[2,2])
    F1  = 100*(2/((1/(SN/100)) + (1/(PPV/100))))

    iresultsT = rbind(iresultsT,c(study,type,lBs[i],SN,SP,PPV,NPV,ACC,F1))
  } # end of i loop
} # end of study loop

iresultsT = data.frame(iresultsT)
colnames(iresultsT) = c("Study", "Type", "thresh", "SN", "SP", "PPV", "NPV", "ACC", "F1")
iresultsT$SN   = as.numeric(as.character(iresultsT$SN))
iresultsT$SP   = as.numeric(as.character(iresultsT$SP))
iresultsT$Type = factor(iresultsT$Type, levels=names(indexes))

write.table(iresultsT, file=paste(analysisdir, "snsp-results.txt", sep=""), col.names=TRUE, row.names=FALSE, sep="\t")

# get study specific AUCplots
aucRes = c()
for (type in levels(iresultsT$Type)){
  aucRow = c(type)
  for (study in levels(iresultsT$Study)){
    Set1   = iresultsT[iresultsT$Type==type & iresultsT$Study==study, ]
    aucEst = trapz(Set1$SP,Set1$SN)/100
    aucRow = c(aucRow, aucEst)
  }
  aucRes = rbind(aucRes,aucRow)
}
colnames(aucRes) = c("Type", levels(iresultsT$Study))

write.table(aucRes, file=paste(analysisdir, "auc-results.txt", sep=""), col.names=TRUE, row.names=FALSE, sep="\t")

colscheme6 = c("#336699",
"#F39C12",
"#F1C40F",
"#27AE60",
"#8E44AD",
"#7F8C8D")

iresultsT$Type = factor(iresultsT$Type, levels=c(
"re-analysis-r" ,   "re-analysis-nr"  , "re-analysis-agg",
"lit-r"       ,     "lit-nr"   ,        "lit-agg",
"lit-melanoma-r"  , "lit-melanoma-nr",  "lit-melanoma-agg"))

outfile1 = paste(pdfdir, "figure.pdf", sep="")
p1<-
ggplot(iresultsT[!grepl("(\\.no\\.|p05)", iresultsT$Type),], aes(x=SP, y=SN), group=Study) +
geom_path(aes(group=Study, color=Study), size=0.5) +
theme_bw() +
scale_colour_manual(values=colscheme6) +
geom_abline(intercept=100, slope=1, color="black", size=0.5, linetype="dashed") +
theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=9, colour="black"),
      axis.text.y = element_text(colour="black", size=9),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x     = element_text(size=7)) +
xlab("Specificity") +
ylab("Sensitivity") +
theme(aspect.ratio = 1) +
scale_x_reverse() +
facet_wrap(~Type, ncol=3)
ggsave(outfile1, plot=p1, width=6, height=6)
