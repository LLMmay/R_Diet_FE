####################################################
library(ggplot2)
data = read.csv("Genus-F.csv",header=T,sep=",")
data$Tissue <- factor(data$Tissue,levels = c("A1","A2","A3","A4","A5","A6","A7","A8","A9","B1","B2","B3","B4","B5","B6","B7","B8","B9","C1","C2"),
                      labels = c("Prevotella","Ruminococcus_E","Butyrivibrio","Limivicinus","UBA1177","UBA1213","UBA1179","Fibrobacter","Treponema_D","Cryptobacteroides","RUG592","Ruminococcus","Hominimerdicola","UBA1711","UBA2450","Pseudobutyrivibrio","RUG14515","Aristaeella","UBA1407","Sodaliphilus"))
a=ggplot(data,aes(x=Tissue,y=Number,fill=KOs))+
  geom_bar(stat = "identity", width = 1,colour="white",size=0.25)+
  coord_polar(theta = "x",start=0) +
  theme_light()+
  scale_fill_manual(values=c("#946769","#4a6d8c"))+   
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey80",size=.25),
        axis.text.y = element_text(size = 11,colour="black"),
        axis.line.y = element_line(size=0.25),
        axis.text.x=element_text(size = 11,colour="black"))
data = read.csv("Genus-G.csv",header=T,sep=",")
data$Tissue <- factor(data$Tissue,levels = c("A1","A2","A3","A4","A5","A6","A7","A8","A9","B1","B2","B3","B4","B5","B6","B7","B8","B9","C1","C2"),
                      labels = c("Prevotella","Selenomonas_A","Cryptobacteroides","Ruminococcus_E","Limivicinus","UBA1067","UBA1179","RUG099","RGIG6816","UBA1248","Faecousia","UBA636","Treponema_D","SIG220","UBA3774","Bact-11","UBA11452","CADBNZ01","UBA1177","Eubacterium_Q"))
b=ggplot(data,aes(x=Tissue,y=Number,fill=KOs))+
  geom_bar(stat = "identity", width = 1,colour="white",size=0.25)+
  coord_polar(theta = "x",start=0) +
  theme_light()+
  scale_fill_manual(values=c("#946769","#4a6d8c"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey80",size=.25),
        axis.text.y = element_text(size = 11,colour="black"),
        axis.line.y = element_line(size=0.25),
        axis.text.x=element_text(size = 11,colour="black"))
library(ggpubr)
library(egg)
ggarrange(a,b,ncol=2)
