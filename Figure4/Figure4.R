library(ggplot2)

library(patchwork)

library(RColorBrewer)

#Plot pathogen-specific effects of rs183868412 genotype

cols <- brewer.pal(11, "RdYlBu")
cols2 <- brewer.pal(8, "Dark2")
cols3 <- brewer.pal(8, "Set2")

label <- rep(c("Aci (n=118)", "BHS (n=130)", "E. coli (n=141)", "Hib (n=113)", "NTS (n=159)", "Other (n=242)","Pneumo (n=391)", "S. aureus (n=151)"),3)
mean  <- c(0.1294529, 0.97544283, 0.94372723, 0.02255999, 0.99779713, 0.80957666, 0.80415958, 0.53290879)
se  <- c(0.4374372, 0.3009951, 0.2873546, 0.4725781, 0.2694641, 0.2438325, 0.2045006, 0.3296354)
lower <- mean-1.96*se
upper <- mean+1.96*se
  
df <- data.frame(label, mean, lower, upper)
df$facet <- "rs183868412:T"

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:8)])
path.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(10,10,10,10)],cols3[8],cols[c(10,10)],cols3[8])) + scale_colour_manual(values = c(cols[c(10,10,10,10)],cols3[8],cols[c(10,10)],cols3[8])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("OR") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(-0.693147181, 0, 0.693147181, 1.386294361),
                     labels=c("0.5",	"1.0",	"2.0",	"4.0"), limits = c(-1,2)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15),plot.title=element_text(size=20))+
  facet_wrap(~facet, ncol=1)+
  ggtitle("Pathogens") +
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold.italic"))

#Plot stratum-specific effects of rs183868412 genotype on neonatal on non-neonatal disease

label <- rep(c("1st 29 days (n=195)", ">28 days (n=1,245)"),3)
mean  <- c(0.9121288, 0.5966586)
se  <- c(0.2547673, 0.1514635)
lower <- mean-1.96*se
upper <- mean+1.96*se

df <- data.frame(label, mean, lower, upper)
df$facet <- "rs183868412:T"

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:2)])
neo.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  ylim(-1,2) +
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(10,10)])) + scale_colour_manual(values = c(cols[c(10,10)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("OR") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(-0.693147181, 0, 0.693147181, 1.386294361),
                     labels=c("0.5",	"1.0",	"2.0",	"4.0"), limits=c(-1,2)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15),plot.title=element_text(size=20))+
  facet_wrap(~facet, ncol=1)+
  ggtitle("Neonatal sepsis") +
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold.italic"))


#Plot stratum-specific effects of rs183868412 genotype on neonatal on non-neonatal disease

label <- rep(c("Bact. no Pf (n=1,236)", "Bact. with Pf (n=204)"),3)
mean  <- c(0.5811352, 0.9540062)
se  <- c(0.1525740, 0.2424632)
lower <- mean-1.96*se
upper <- mean+1.96*se

df <- data.frame(label, mean, lower, upper)
df$facet <- "rs183868412:T"

levels(df$label)

cols3 <- brewer.pal(8, "Set2")


# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label)[c(1:2)])
para.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[c(10,10)])) + scale_colour_manual(values = c(cols[c(10,10)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("OR") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(-0.693147181, 0, 0.693147181, 1.386294361),
                     labels=c("0.5",	"1.0",	"2.0",	"4.0"), limits=c(-1,2)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15),plot.title=element_text(size=20))+
  facet_wrap(~facet, ncol=1)+
  ggtitle("Parasitaemia") +
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold.italic"))

#Plot stratum-specific effects of rs183868412 genotype on risk of bacteraemia over study period

cols <- brewer.pal(11, "RdYlBu")
cols2 <- brewer.pal(8, "Dark2")


label <- c("95", "96", "97", "98", "98/99/00", "00","01/02","02","03/04","04","05/06","06","07/08","08", "09/10", "10","11" , "12", "13","14", "15")
mean  <- c(NA,NA,NA,NA,0.826615, NA, 0.796132,NA,0.546421,NA,0.510647,NA,0.567562,NA,1.03668,NA,NA,NA,NA,NA,NA)
se  <- c(NA,NA,NA,NA,0.174712, NA,0.203987,NA, 0.200618, NA, 0.242699,NA, 0.268059,NA, 0.509352,NA,NA,NA,NA,NA,NA)
lower <- mean-1.96*se
upper <- mean+1.96*se

df <- data.frame(label, mean, lower, upper)
df$facet <- "Bacteraemia risk & rs183868412:T"

levels(df$label)

cols3 <- brewer.pal(8, "Set2")


# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=c("95", "96", "97", "98", "98/99/00", "00","01/02","02","03/04","04","05/06","06","07/08","08", "09/10", "10","11" , "12", "13","14", "15"))
yr.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  scale_x_discrete(breaks=c("95","00","05/06","10"),
                                 labels=c("1995", "2000", "2005", "2010"),expand = expand_scale(add = 1))+
  geom_pointrange(fatten=4, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = cols[rep(10,21)]) + scale_colour_manual(values = c(cols[rep(10,21)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  ylab("OR") +
  scale_y_continuous(breaks=c(-0.693147181, 0, 0.693147181, 1.386294361, 2.079441542),
                     labels=c("0.5",	"1.0",	"2.0",	"4.0", "8.0"), limits=c(-1,2.5)) +
  theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15),panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())+
  facet_wrap(~facet, ncol=1)+
  theme(strip.background = element_rect(color="black", fill=cols3[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold.italic"))

pp = read.table("/Users/jamesgilchrist/Documents/g6pd/parasite_prev/para_prev.txt", header = TRUE)

#Plot falling parasite prevalence among trauma admission over the study period

para_prev<-ggplot(data=pp, aes(x=year, y=malinc)) + geom_line(lwd=1.5) +
geom_ribbon(aes(ymin=low, ymax=high), linetype=2, alpha=0.3, fill = "grey70", color = "red", 
            outline.type = "both") +
  theme_bw() + xlab("Year") + ylab("Prevalence (%)") + theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank(), 
                                                             axis.title=element_text(size=15), 
                                                             axis.text=element_text(size=12),
                                                             plot.title=element_text(size=18))

Figure4 <- (path.fp|(neo.fp/para.fp))/(yr.fp/para_prev) + plot_annotation(tag_levels = 'A')

