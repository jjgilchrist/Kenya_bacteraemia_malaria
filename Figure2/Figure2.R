library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(ggplotify)
library(patchwork)

#thinned example summary pvalues 1:100 SNPs

disc.add <- read.table("example.pval.gz", header = F)
colnames(disc.add) <- c("chr", "bp", "rsid", "p")
str(disc.add)

disc.add$chr <- factor(disc.add$chr, levels = c(1:22))
levels(disc.add$chr) <- c(1:22)

#highlight chr2 hit
disc.add$label <- NA
disc.add$label[which(disc.add$rsid=="rs183868412:T")] <- "BIRC6"

 
disc.add$anno <- NA
disc.add$anno[which(disc.add$rsid=="rs183868412:T")] <- 1

#data for manahattan
don.1 <- disc.add %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  left_join(disc.add, ., by=c("chr"="chr")) %>%
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = don.1 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

#Manhattan plot
Figure2 <- ggplot(don.1, aes(x=BPcum, y=-log10(p))) +
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  geom_text_repel( data=subset(don.1, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = don.1, linetype="dashed", color = "red") +
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","2","4","6","8","10"), breaks = c(0,2,4,6,8,10), expand = c(0, 0), limits= c(0,12)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

#data for QQ plot
df <- data.frame(observed = -log10(sort(na.omit(disc.add$p))),
    expected = -log10(ppoints(length(na.omit(disc.add$p)))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  chisq <- qchisq(1-disc.add$p,1)
  median(chisq)/qchisq(0.5,1)
  
  lambda <- paste("lambda == ", round(median(chisq)/qchisq(0.5,1),3))

#QQ plot
Figure2.supp4 <- ggplot(df) +
  geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
  geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
  xlab(log10Pe) +
  ylab(log10Po) +
  annotate("text", x = 1, y = 9, label = lambda, parse = TRUE, size = 10) +
  theme_bw() + ylim(NA, 10) + 
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))


#Genotyping concordance data
concord <- read.table("genotyping.concord.txt", header = TRUE)


qplot(concord$pairwise,
      geom="histogram",
      binwidth = 0.001,  
      main = "genotype concordance", 
      xlab = "pairwise concordance",
      ylab = "count",
      xlim = c(0.8,1.01))

median(concord$pairwise)

#Genotyping concordance plot

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

Figure2.supp1 <- ggplot(data=concord, aes(pairwise)) + 
  geom_histogram(breaks=seq(0.8, 1, by = 0.001), col = cols3[2]) + 
  labs(x="pairwise concordance", y="count") + 
  xlim(c(0.8,1.01)) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20))

#Principal components data
total <- read.table("PCs.txt", header = T)

#Principal components: genotyping platform plots

p.pc12 = ggplot(total, aes(PC1, PC2)) + 
  geom_point(data=subset(total, PLATFORM=="AFFY"), color=cols[4], size=4, alpha = 0.15) +
  geom_point(data=subset(total, PLATFORM=="OMNI"), color=cols[5], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC1") +
  xlab("PC2") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

p.pc34 = ggplot(total, aes(PC3, PC4)) + 
  geom_point(data=subset(total, PLATFORM=="AFFY"), color=cols[4], size=4, alpha = 0.15) +
  geom_point(data=subset(total, PLATFORM=="OMNI"), color=cols[5], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC3") +
  xlab("PC4") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

p.pc56 = ggplot(total, aes(PC5, PC6)) + 
  geom_point(data=subset(total, PLATFORM=="AFFY"), color=cols[4], size=4, alpha = 0.15) +
  geom_point(data=subset(total, PLATFORM=="OMNI"), color=cols[5], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC5") +
  xlab("PC6") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

Figure2.supp2 <- (p.pc12| p.pc34| p.pc56)



#Principal components: self-reported ethnicity plots

p.pc12 = ggplot(total, aes(PC1, PC2)) + 
  geom_point(data=subset(total, eth==4), color=cols[9], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==1), color=cols[1], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==2), color=cols[2], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==3), color=cols[3], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC1") +
  xlab("PC2") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

p.pc34 = ggplot(total, aes(PC3, PC4)) + 
  geom_point(data=subset(total, eth==4), color=cols[9], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==1), color=cols[1], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==2), color=cols[2], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==3), color=cols[3], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC3") +
  xlab("PC4") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))

p.pc56 = ggplot(total, aes(PC5, PC6)) + 
  geom_point(data=subset(total, eth==4), color=cols[9], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==1), color=cols[1], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==2), color=cols[2], size=4, alpha = 0.15) +
  geom_point(data=subset(total, eth==3), color=cols[3], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC5") +
  xlab("PC6") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))


Figure2.supp3 <- (p.pc12| p.pc34| p.pc56)
