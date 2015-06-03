##################### Bayescan ########################3
##### Script R

### Script bayescan

###################################### REMOVE BALANCING SELECTION ###############################################

#setwd("/Volumes/thierry_mac/populations_esturgeon/adu_allpop_m10_r70/MAF5/bayescan/")

bayescan_file=read.table("input_file.temp",dec=".")
SNPb=read.table("../../01-loci/pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17-keep-locis.txt",header=T)
bayescan=cbind(SNPb, bayescan_file)
colnames(bayescan)=c("SNP","POST_PROB","LOG10_PO","Q_VALUE","ALPHA","FST")

# POST_PROB = 1 & Q_VALUE = 0 == 0.0001 (see comment below on LOG10_PO)
attach(bayescan)
class(bayescan$Q_VALUE)
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE)
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001

bayescan$POST_PROB <- (round(bayescan$POST_PROB, 4))
bayescan$LOG10_PO <- (round(bayescan$LOG10_PO, 4))
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4))
bayescan$ALPHA <- (round(bayescan$ALPHA, 4))
bayescan$FST <- (round(bayescan$FST, 6))
rm(bayescan_file)

# ADD a column for selection grouping
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing"))
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

positive <- bayescan[bayescan$SELECTION=="diversifying",]
neutral <- bayescan[bayescan$SELECTION=="neutral",]
balancing <- bayescan[bayescan$SELECTION=="balancing",]
xtabs(data=bayescan, ~SELECTION)

#Transformation Log de la valeur Q pour le graph
bayescan$LOG10_Q <- log10(bayescan$Q_VALUE)

#Quantile of FST
FST_quantile<-quantile(bayescan$FST,c(0,0.2,0.4,0.6,0.8,1.0),na.rm="true")
bayescan$FST_GROUP <- cut(bayescan$FST, FST_quantile,labels=c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))

# Grouping des groupes LOG10_PO
bayescan$PO_GROUP <- ifelse(bayescan$LOG10_PO>2,"decisive",ifelse(bayescan$LOG10_PO>1.5,"very strong",ifelse(bayescan$LOG10_PO>1,"strong",ifelse(bayescan$LOG10_PO>0.5,"substantial","no evidence"))))
bayescan$PO_GROUP<- factor(bayescan$PO_GROUP)
levels(bayescan$PO_GROUP)
class(bayescan$PO_GROUP)
bayescan[["PO_GROUP"]] <- factor(bayescan[["PO_GROUP"]],levels=c("no evidence","substantial","strong","very strong","decisive"),ordered=T)
levels(bayescan$PO_GROUP)

# To have an idea in numbers...
xtabs(data=bayescan, ~PO_GROUP)

##TODO change name to figure
#################################### SCATTER PLOT OF BAYESCAN RESULTS ###########################################
library(ggplot2)

graph_title="Bayescan plot for 14 276 SNP"
x_title="Log10(Q_VALUE)"
y_title="Fst"
posi=length(positive$SELECTION)
neut=length(neutral$SELECTION)
bal=length(balancing$SELECTION)

graph_1<-ggplot(bayescan,aes(x=LOG10_Q,y=FST))
graph_1+geom_point(aes(colour=bayescan$PO_GROUP,shape=FST_GROUP))+
  scale_shape_manual(name="FST quantile\ngroup",values=c(5,2,3,4,1))+
  scale_colour_manual(name="Model choice",values=c("darkred","yellow","orange","green","forestgreen"))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  geom_vline(xintercept=c(log10(0.05)),color="darkgreen")+
  annotate("text", x= -2.7, y=0.125,label=posi,colour="green")+
  annotate("text", x= -0.7, y=0.15,label=bal,colour="red")+
  annotate("text", x= -0.7, y=0.07,label=neut,colour="black")+
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"))
ggsave("bayescan_anguilla-ecotype.png")
dev.off()

#### Select best loci####

write.table(bayescan, "recode-bayescan.txt", row.names=F, quote=F)
## Qvalue<=0.6
best1=subset(bayescan,subset=(bayescan$Q_VALUE<=0.6))
length(best1$SNP)
write.table(best1$SNP, "Q6.txt",row.names=F, quote=F)

## Qvalue<=0.8

#best2=subset(bayescan,subset=(bayescan$Q_VALUE<=0.8))
#length(best2$SNP)
#write.table(best2$SNP, "pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.Q8.txt" ,row.names=F, quote=F)
