library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(patchwork)
library(nnet)

for_plot <- read.table("birc6_for_RA_b38.txt", header = T)

birc6.region <- subset(for_plot, bp>31277134 & bp<33274933)

birc6.region$bin_r2 <- 1
birc6.region$bin_r2[which(birc6.region$r2>0.1 & birc6.region$r2 <= 0.3)] <- 2
birc6.region$bin_r2[which(birc6.region$r2>0.3 & birc6.region$r2 <= 0.5)] <- 3
birc6.region$bin_r2[which(birc6.region$r2>0.5 & birc6.region$r2 <= 0.7)] <- 4
birc6.region$bin_r2[which(birc6.region$r2>0.7)] <- 5

birc6.region$bin_dp <- 1
birc6.region$bin_dp[which(birc6.region$dp>0.6 & birc6.region$dp <= 0.7)] <- 2
birc6.region$bin_dp[which(birc6.region$dp>0.7 & birc6.region$dp <= 0.8)] <- 3
birc6.region$bin_dp[which(birc6.region$dp>0.8 & birc6.region$dp <= 0.9)] <- 4
birc6.region$bin_dp[which(birc6.region$dp>0.9)] <- 5

birc6.region$annotate <- 0
birc6.region$annotate[c(which(birc6.region$rsid=="rs183868412:T"))] <- 1
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(8,"Set2")

#BIRC6 - regional assocaition plot (discovery data)

birc6_total <- ggplot(birc6.region, aes(x=bp, y=-log10(p))) + 
  xlim(31874931, 32874933) +
  geom_point(data=subset(birc6.region, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,2,4,6,8,10)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(birc6.region, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_y = 1) +
  annotate("text", x = (32180000-225069), y = 11, label = quote(r^2), hjust = 0.5) +
  annotate("point", x = (32120000-225069), y = 10.3, size = 3, colour = cols2[1]) +
  annotate("text", x = (32180000-225069), y = 10.3, label = c("0.7-1.0"), hjust = 0.5) +
  annotate("point", x = (32120000-225069), y = 9.6, size = 3, colour = cols2[2]) +
  annotate("text", x = (32180000-225069), y = 9.6, label = c("0.5-0.7"), hjust = 0.5) +
  annotate("point", x = (32120000-225069), y = 8.9, size = 3, colour = cols2[3]) +
  annotate("text", x = (32180000-225069), y = 8.9, label = c("0.3-0.5"), hjust = 0.5) +
  annotate("point", x = (32120000-225069), y = 8.2, size = 3, colour = cols2[5]) +
  annotate("text", x = (32180000-225069), y = 8.2, label = c("0.1-0.3"), hjust = 0.5) +
  annotate("point", x = (32120000-225069), y = 7.5, size = 3, colour = cols[8]) +
  annotate("text", x = (32180000-225069), y = 7.5, label = c("<0.1"), hjust = 0.5) +
  annotate("text", x = (33100000-225069), y = 11, label = c("Discovery"), hjust = 1, size=10)

#BIRC6 - regional assocaition plot (replication data)

birc6_rep <- ggplot(birc6.region, aes(x=bp, y=-log10(p_rep))) + 
  xlim(31874931, 32874933) +
  geom_point(data=subset(birc6.region, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==4), color=cols2[2], size=3) + 
  geom_point(data=subset(birc6.region, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log10(p)", breaks=c(0,2,4), limits=c(NA,4)) + 
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(birc6.region, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_y = 0.5) +
  annotate("text", x = 32874933, y = 3.8, label = c("Replication"), hjust = 1, size=10)


#plot genes
genes <- read.table("genes_b38.txt", header = T)

plot.range <- c(31874931, 32874933)
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(31874931, 32874933) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes[seq(9,15,2),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(10,16,2),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(9,15,2),], aes(x=mid, label=external_gene_name), size=4, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(10,16,2),], aes(x=mid, label=external_gene_name), size=4, col = c("black"),
                   nudge_y =0.05, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#plot recombination rate

recomb <- read.table("genetic_map_chr2_b38.txt", header = F)
colnames(recomb) <- c("position", "rate", "map")
recomb_birc6 <- subset(recomb, position>31874931 & position<32874933)

recomb_rate <- ggplot(recomb_birc6, aes(x=position, y=rate)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 2") +
  scale_x_continuous(breaks=c(32000000,32200000,32400000,32600000,32800000),
                     labels=c("32.0Mb", "32.2Mb", "32.4Mb", "32.6Mb", "32.8Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))



cols <- brewer.pal(11,"RdYlBu")

data <- read.table("prob_plot.txt", header = T)

#Plot distribution of 'true' severe malaria probabilities

data$prob_strata_half <- "SM: P(SM|Data)<0.5"
data$prob_strata_half[which(data$p_bact<0.5)] <- "SM: P(SM|Data)>0.5"
data$prob_strata_half[which(data$prob_strata=="BACT")] <- "BACT"
data$prob_strata_half[which(data$prob_strata=="CONTROL")] <- "CONTROL"
data$prob_strata_half <- factor(data$prob_strata_half)
data$prob_strata_half <- relevel(data$prob_strata_half, ref = "CONTROL")

data1 <- data[-c(which(data$prob_strata_half=="CONTROL"), which(data$prob_strata_half=="BACT")),]

prob_SM_hist<-ggplot(data1, aes(x=(1-p_bact), fill=factor(prob_strata_half), colour=factor(prob_strata_half))) + 
  geom_histogram(binwidth=0.012, alpha=0.72) +
  scale_fill_manual(values = cols[c(4,2)]) +
  scale_colour_manual(values = cols[c(4,2)]) +
  theme_bw() + xlab("P(Severe Malaria | Data)") + theme(legend.position = "none") +
  geom_hline(yintercept=0, color=cols[8]) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

cols <- brewer.pal(8, "Set2")
cols2 <- brewer.pal(8, "Dark2")

#Plot effect of rs86868412 genotype across case strata

label <- rep(c("SM - P.high", "SM - P.low", "Bacteraemia - disc", "Bacteraemia - rep"),3)
mean  <- c(
  0.07180831, 
  0.8637545, 
  0.7536469, 
  1.06025)
se  <- c(
  0.3202724, 
  0.3192534, 
  0.1447982, 
  0.313471)
lower <- mean-1.96*se
upper <- mean+1.96*se

df <- data.frame(label, mean, lower, upper)
df$facet <- "rs183868412:T"

levels(df$label)

cols3 <- brewer.pal(8, "Set2")
cols <- brewer.pal(11,"RdYlBu")

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:4)])
birc6.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = cols[c(10,10,4,2)]) + scale_colour_manual(values = cols[c(10,10,4,2)]) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("OR") + 
  scale_y_continuous(breaks=c(-0.693147181, 0, 0.693147181, 1.386294361),
                     labels=c("0.5",	"1.0",	"2.0",	"4.0"), limits=c(-1,2)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15))+
  facet_wrap(~facet, ncol=1)+
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"))

Figure3 <- ((birc6_total/birc6_rep/genes.plot/recomb_rate)/(prob_SM_hist|birc6.fp)) + plot_layout(heights = c(3, 1.5, 1, 1,3))


#Figure 3 Supplement 1: plot principal components of genotyping data vs. self-reported ethnicity
cols <- brewer.pal(9,"Set1")
total <- read.table("pop_plots_rep.txt", header = T)

p.pc12 = ggplot(total, aes(PC1, PC2)) + 
  geom_point(data=subset(total, ETHNICITY==4), color=cols[9], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==1), color=cols[1], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==2), color=cols[2], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==3), color=cols[3], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC1") +
  xlab("PC2") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
p.pc12

p.pc34 = ggplot(total, aes(PC3, PC4)) + 
  geom_point(data=subset(total, ETHNICITY==4), color=cols[9], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==1), color=cols[1], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==2), color=cols[2], size=4, alpha = 0.15) +
  geom_point(data=subset(total, ETHNICITY==3), color=cols[3], size=4, alpha = 0.15) +
  theme_bw() +
  ylab("PC3") +
  xlab("PC4") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=25))
p.pc34


Figure3.Supp1 <- (p.pc12| p.pc34)



