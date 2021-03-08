
  rm(list=ls())
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)

  setwd(paste0(path.expand("~"),"/Desktop/shaikh-sears-clin-cancer-research-2021/code/"))

  # ------------------
  analysisdir = "../analysis/figureS1/"
  unlink(analysisdir, recursive=TRUE)
  dir.create(analysisdir)
  # ------------------

  A <- read.table("../data/final/pre-S1/results.txt", header=TRUE, sep='\t', check.names=FALSE)
  A$Study = factor(A$Study, levels=c("Gopalakrishnan","Matson","Frankel"))
  A = A[sample(1:nrow(A)),]

  Matson         = A[A$Study=="Matson",]
  Gopalakrishnan = A[A$Study=="Gopalakrishnan",]
  Frankel        = A[A$Study=="Frankel",]
  mp <- cor.test(Matson$abWGS, Matson$ab16SrRNA,method="pearson")
  ms <- cor.test(Matson$abWGS, Matson$ab16SrRNA,method="spearman")
  gp <- cor.test(Gopalakrishnan$abWGS, Gopalakrishnan$ab16SrRNA,method="pearson")
  gs <- cor.test(Gopalakrishnan$abWGS, Gopalakrishnan$ab16SrRNA,method="spearman")
  fp <- cor.test(Frankel$abWGS, Frankel$ab16SrRNA,method="pearson")
  fs <- cor.test(Frankel$abWGS, Frankel$ab16SrRNA,method="spearman")

  result = c(mp$estimate, ms$estimate, gp$estimate, gs$estimate, fp$estimate, fs$estimate, mp$p.value, ms$p.value, gp$p.value, gs$p.value, fp$p.value, fs$p.value)

  colscheme3 = c("#E86850",
                 "#447294",
                 "#ddbc4c",
                 "#ccc9c9")

  outfile1 = paste0(analysisdir,"plot1.pdf")
  p1 <- ggplot(A, aes(x=abWGS, y=ab16SrRNA)) +
        geom_point(aes(colour = Study), size=2.0, alpha=0.85) +
        scale_colour_manual(values=colscheme3, drop=FALSE) +
        theme_bw() + theme(aspect.ratio = 1) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=12, colour="black"),
              axis.text.y = element_text(colour="black", size=12)) +
        xlab("WGS abundance") +
        ylab("16S rRNA abundance") +
        scale_x_log10(breaks=c(0.01, 0.1, 1, 2, 4, 10, 25, 50, 100), minor_breaks=NULL) +
        scale_y_log10(breaks=c(0.01, 0.1, 1, 2, 4, 10, 25, 50, 100), minor_breaks=NULL) +
        theme(legend.position="top") +
        geom_abline(slope=1, intercept=0, color="black", linetype = "dashed", alpha=0.5)
  ggsave(outfile1,p1)

  outfile1 = paste0(analysisdir,"plot2.pdf")
  p1 <- ggplot(A, aes(x=abWGS, y=ab16SrRNA)) +
  geom_point(aes(colour = Study), size=1.0, alpha=0.75) +
  scale_colour_manual(values=colscheme3, drop=FALSE) +
  theme_bw() + theme(aspect.ratio = 0.8) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6, colour="black"), axis.text.y = element_text(colour="black", size=6), strip.text.x = element_text(colour="black", size=7)) +
  xlab("WGS abundance") +
  ylab("16S rRNA abundance") +
  scale_x_log10(breaks=c(0.01, 0.1, 1, 2, 4, 10, 25, 50, 100), minor_breaks=NULL) +
  scale_y_log10(breaks=c(0.01, 0.1, 1, 2, 4, 10, 25, 50, 100), minor_breaks=NULL) +
  facet_wrap(~Species) +
  theme(legend.position="top") +
  geom_abline(slope=1, intercept=0, color="black", linetype = "dashed", alpha=0.5)
  ggsave(outfile1,p1)
