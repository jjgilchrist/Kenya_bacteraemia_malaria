library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)


cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(11,"RdYlBu")

ps <- read.table("Fig1_data.txt", header = T)

#Plot PfHRP2 concentration vs platelet count - colour by probability 'true' severe malaria, i.e. Model 1

prob_hrp2 = ggplot(ps, aes(log10(platelet), log10(hrp2))) + 
  geom_point(aes(colour = factor(p2_fifths)), size = 2) + scale_colour_manual(values =cols3[c(5,4,3,2,1)], na.translate=FALSE, name = "P(SM | Data)", labels = c("≤0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8")) +
  theme_bw() +
  scale_y_continuous(breaks=c(0,1,2,3,4,5),
                   labels=c(1,10,expression(10^2), expression(10^3), expression(10^4), expression(10^5))) +
  scale_x_continuous(breaks=c(log10(10),log10(20), log10(50), log10(200), log10(1000)),
                     labels=c(10,20,50,200,1000)) +
  labs(x = expression(paste("Platelet count (x", 10^9,"/L)")),
        y = "PfHRP2 (ng/ml)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        plot.title=element_text(size=50, face = "bold"))

#Plot white cell count vs platelet count - colour by probability 'true' severe malaria, i.e. Model 2

prob_wcc = ggplot(ps, aes(log10(platelet), log10(wbc))) + 
  geom_point(aes(colour = factor(p1_fifths)), size = 2) + scale_colour_manual(values =cols3[c(5,4,3,2,1)], na.translate=FALSE, name = "P(SM | Data)", labels = c("≤0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8")) +
  theme_bw() +
  scale_y_continuous(breaks=c(log10(1),log10(2),log10(5), log10(10), log10(20), log10(50), log10(100)),
                     labels=c(1,2,5,10, 20,50,100)) +
  scale_x_continuous(breaks=c(log10(10),log10(20), log10(50), log10(200), log10(1000)),
                     labels=c(10,20,50,200,1000)) +
  labs(x = expression(paste("Platelet count (x", 10^9,"/L)")),
       y = expression(paste("White cell count (x", 10^9,"/L)"))) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        plot.title=element_text(size=50, face = "bold"))

cols <- brewer.pal(11,"RdYlBu")

#Distribution of 'true' malaria probabilities: Model 1

hist_hrp2<-ggplot(ps, aes(x=P_SM2, fill=factor(p2_fifths), colour=factor(p2_fifths))) + 
  geom_histogram(binwidth=0.012, alpha=0.72) +
  scale_fill_manual(values = cols[c(5,4,3,2,1)]) +
  scale_colour_manual(values = cols[c(5,4,3,2,1)]) +
  theme_bw() + xlab("P(SM | Data)") + theme(legend.position = "none") +
  geom_hline(yintercept=0, color=cols[8]) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

#Distribution of 'true' malaria probabilities: Model 2

hist_wcc<-ggplot(ps, aes(x=P_SM1, fill=factor(p1_fifths), colour=factor(p1_fifths))) + 
  geom_histogram(binwidth=0.012, alpha=0.72) +
  scale_fill_manual(values = cols[c(5,4,3,2,1)]) +
  scale_colour_manual(values = cols[c(5,4,3,2,1)]) +
  theme_bw() + xlab("P(SM | Data)") + theme(legend.position = "none") +
  geom_hline(yintercept=0, color=cols[8]) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(11,"RdYlBu")

#PfHRP2 concentration as a predictor of bacteraemia

p_labels = data.frame(expt = c("bloodculture_pos"), label = c("italic(P)==9.62 %*% 10^-6"))
hrp2.box = ggplot(ps, aes(x=factor(bloodculture_pos), y=log10(hrp2))) +
geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
scale_x_discrete(breaks=c(0,1), labels = c("No", "Yes")) + aes(fill = factor(bloodculture_pos), col=factor(bloodculture_pos)) + scale_fill_manual(values = cols[c(8,3)]) + scale_colour_manual(values = cols[c(8,3)]) +
ylab("PfHRP2 (ng/ml)") + xlab("Bacteraemia") +
scale_y_continuous(breaks=c(0,1,2,3,4,5),
                     labels=c(1,10,expression(10^2), expression(10^3), expression(10^4), expression(10^5))) +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + 
  geom_text(x = 1.5, y = Inf, vjust=1.5, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)


#Platelet count as a predictor of bacteraemia

p_labels = data.frame(expt = c("bloodculture_pos"), label = c("italic(P)==0.0144"))
plt.box = ggplot(ps, aes(x=factor(bloodculture_pos), y=log10(platelet))) +
  geom_dotplot(binaxis="y", binwidth=0.01, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  scale_x_discrete(breaks=c(0,1), labels = c("No", "Yes")) + aes(fill = factor(bloodculture_pos), col=factor(bloodculture_pos)) + scale_fill_manual(values = cols[c(8,3)]) + scale_colour_manual(values = cols[c(8,3)]) +
  ylab(expression(paste("Platelet count (x", 10^9,"/L)"))) + xlab("Bacteraemia") +
  scale_y_continuous(breaks=c(log10(10),log10(20), log10(50), log10(200), log10(1000)),
                     labels=c(10,20,50,200,1000)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) + 
  geom_text(x = 1.5, y = Inf, vjust=2, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)

#White cell count as a predictor of bacteraemia

p_labels = data.frame(expt = c("bloodculture_pos"), label = c("NS"))
wcc.box = ggplot(ps, aes(x=factor(bloodculture_pos), y=log10(wbc))) +
  geom_dotplot(binaxis="y", binwidth=0.01, stackdir="center", alpha = 0.75) + geom_boxplot(alpha = 0.5) +
  scale_x_discrete(breaks=c(0,1), labels = c("No", "Yes")) + aes(fill = factor(bloodculture_pos), col=factor(bloodculture_pos)) + scale_fill_manual(values = cols[c(8,3)]) + scale_colour_manual(values = cols[c(8,3)]) +
  ylab(expression(paste("White cell count (x", 10^9,"/L)"))) + xlab("Bacteraemia") +
  scale_y_continuous(breaks=c(log10(1),log10(2),log10(5), log10(10), log10(20), log10(50), log10(100)),
                     labels=c(1,2,5,10, 20,50,100)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) + 
  geom_text(x = 1.5, y = Inf, vjust=2, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)


Figure1 <- ((prob_hrp2|hist_hrp2)/(plt.box|hrp2.box)) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

Figure1.supp1 <- ((prob_wcc|hist_wcc)/(plt.box|wcc.box))+ plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

