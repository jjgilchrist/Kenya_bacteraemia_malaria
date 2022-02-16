library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(coloc)


#Colocalisation of BIRC6 exon eQTL with invasive bacterial disease association at rs183868412
d1 <- read.table("exon_birc6.txt", header = T)
ld <- as.matrix(read.table("ld_for_coloc.txt", row.names = 1, header = T))

eqtl.list <- list(beta = d1$BETA, varbeta = (d1$SE)^2, type = "quant", N = 200, MAF = d1$maf, snp = as.character(d1$SNP), LD=ld)
gwas.list <- list(beta = d1$gwas.beta, varbeta = (d1$gwas.se)^2, type = "cc", N = 5400, s = 0.479, snp = as.character(d1$SNP), LD=ld)
E=runsusie(eqtl.list,coverage=0.01)
G=runsusie(gwas.list, coverage=0.01)
susie.res=coloc.susie(E,G)
h4 <- round(max(susie.res$summary[,8]),2)
hit1 <- as.character(susie.res$summary[which(susie.res$summary[,8]==max(susie.res$summary[,8])),2]$hit1)
hit2 <- as.character(susie.res$summary[which(susie.res$summary[,8]==max(susie.res$summary[,8])),3]$hit2)

#plot colocalsation
exon <- d1

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"Spectral")

exon$bin_r2 <- 1
exon$bin_r2[which(exon$r2>0.1 & exon$r2 <= 0.3)] <- 2
exon$bin_r2[which(exon$r2>0.3 & exon$r2 <= 0.5)] <- 3
exon$bin_r2[which(exon$r2>0.5 & exon$r2 <= 0.7)] <- 4
exon$bin_r2[which(exon$r2>0.7)] <- 5

exon.plot <- ggplot(exon, aes(x=-log10(gwas.p), y=-log10(PVAL))) + 
  geom_point(data=subset(exon, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(exon, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(exon, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(exon, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(exon, bin_r2==5), color=cols2[1], size=3) +
  ylim(0,5) +
  xlim(0,10) +
  ylab("BIRC6 exon eQTL -log10(p)") + 
  xlab("Bacterial disease GWAS -log10(p)") + 
  annotate("text", x = 0, y = 5, label = paste0("PP4=",h4), size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), plot.title=element_text(size=20))


#Plot effects of rs184868412 genotype on exon expression in monocytes across stimulation contexts

label <- c("Naive", "IAV", "LPS", "Pam3CSK4", "R848")
mean  <- c(-0.155186,-0.0536451, 0.227218, -0.654215, -0.386583)
se  <- c(0.138402, 0.248336, 0.152519, 0.156258, 0.198497)
lower <- mean-1.96*se
upper <- mean+1.96*se

df <- data.frame(label, mean, lower, upper)
df$facet <- "rs183868412:T"

cols3 <- brewer.pal(8, "Set2")
df$label <- factor(df$label, levels=rev(df$label)[c(3,1,2,4,5)])

eqtl.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols3[c(8,8,3,8,8)])) + scale_colour_manual(values = c(cols3[c(8,8,3,8,8)])) +
  geom_hline(yintercept=0, lty=2) +  
  coord_flip() +  
  ylab("beta") + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15),plot.title=element_text(size=20))+
  facet_wrap(~facet, ncol=1)+
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold.italic"))

Figure5 <- (eqtl.fp|exon.plot) + plot_annotation(tag_levels = 'A')


