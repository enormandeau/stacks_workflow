##################### Run king ################################
bed_file=pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17

# PCA
king -b $bed_file.bed --pca --prefix $bed_file"_pca_"

# MDS
king -b $bed_file.bed --mds --prefix $bed_file"_mds_"

#Allele Frequency Statistics

king -b $bed_file.bed --individual --prefix $bed_file"_ind_"

# relationship

king -b $bed_file.bed --kinship --prefix $bed_file"_relationship_"

################### CODE R ######################################
library(ggplot2)

king=read.table("pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17_pca_pc.ped", header = F)
names(data)

#PC1_2
graph_title="ACP in Anguilla anguilla"
n=nrow(king) # TO ANNOTATE THE SNP NUMBER ON GRAPH
x_title="PC5"
y_title="PC6"
graph_1<-ggplot(king,aes(x=king$V11,y=king$V12))
graph_1+geom_point()
 graph_1+geom_point(aes(colour=king$V1))+
 scale_colour_manual(name="Populations",values=rainbow(16))+
 labs(title=graph_title)+
 labs(x=x_title)+
 labs(y=y_title)

ggsave("anguilla_acp_king_PCA5_6.pdf",width=7,height=5,dpi=600)
dev.off()


#### Eliminer SNP en déséquilibre ####

vcftools --vcf pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf --geno-r2

### R ####

ld=read.table("out.geno.ld", header=T)
dif=abs(ld$POS1-ld$POS2)
lddif=data.frame(ld, dif)
plot(lddif$dif, lddif$R.2)

ld1=lddif[lddif$dif<=90,]
length(ld1$dif)
ld2=ld1[ld1$R.2>=0.2,]
length(ld2$dif)
plot(ld1$dif, ld1$R.2)
write.table(ld2, "locis-desequilibre.txt",row.names=F, quote=F)

### Linux
cut -d " " -f 1,2 locis-desequilibre.txt | sort| uniq| perl -pe "s/\|.+ | /-/g" > locis-desequilibre-uniq.txt

vcftools --vcf pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf --exclude locis-desequilibre-uniq.txt --out ecotype-LD --recode
vcftools --vcf ecotype-LD.recode.vcf --plink-tped --out ecotype-LD
perl -i.bak -pe "s/^A\..-//g|s/-[0-9]+.A/\tA/g" ecotype-LD.tfam 
plink --noweb --tfile ecotype-LD --recode12 --out ecotype-LD
plink --noweb --tfile ecotype-LD --recodeA --out ecotype-LD

#################### Run Adegenet ###########################

library(adegenet)
library(StAMPP)
library(ade4)
library(pegas)
library(ggplot2)

data5=read.PLINK("pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.raw", map.file="pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.map", quiet=FALSE, chunkSize=1000,,parallel=FALSE, n.cores=NULL)


### Statistique générale

na=glNA(data5)
write.table(na, "missing-par-locus.txt")

### frequence de l'allele mineur
myFreq=glMean(data5)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies", breaks=100, main="Distribution of (second) allele frequencies")


#### PCA adegenet par individus
pca5=glPca(data5)

scatter(pca, posi="bottomright")
title("PCA of Anguilla axes 1-2")

# ecrire resultats
col9=data.frame(pca9$scores)
load9=data.frame(data9$loc.names,pca9$loadings)
write.table(load9, "loadings-pca-2pop-ecotype-LD.txt",row.names=F)
write.table(pca9$scores, "scores-pca-2pop-ecotype-LD.txt")

### information ecotype
eco=read.table("../../../16pop-maf0.1/06-admixture_plink/pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.tfam")
eco$V6=factor(eco$V6)

#PC1_2
graph_title="ACP in Anguilla anguilla"
n=nrow(col9) # TO ANNOTATE THE SNP NUMBER ON GRAPH
x_title="PC1"
y_title="PC2"
graph_1<-ggplot(col9,aes(x=col9$PC1 ,y=col9$PC2))
graph_1+geom_point()
graph_1+geom_point(aes(colour=data5$pop, shape=eco$V6))+
scale_colour_manual(name="Population",values=rainbow(16))+
scale_shape_manual(name="Ecotype",values=c(2,3))+
labs(title=graph_title)+
labs(x=x_title)+
labs(y=y_title)


ggsave("anguilla_acp_2popecotype_PCA1_2.pdf",width=10,height=5,dpi=600)
dev.off()



### frequence de l'allele mineur pour les meilleurs SNP

best=load5[load5$Axis1>=0.02 | load5$Axis1<=-0.02 ,]
length(best$Axis1)
bestfreq=myFreq[names(myFreq) %in% best$data5.loc.names]
hist(bestfreq, proba=TRUE, col="gold", xlab="Allele frequencies", breaks=20, main="Distribution of Delta p frequencies")
write.table(best$data5.loc.names, "SNP-best0.02.txt", row.names=F, quote=F)

### changement de freq entre groupe da la PCA
#faire les froupes selon score PCA
score=data.frame(data5$ind.names, col5)
grp1=score[score$PC1>=0,]
length(grp1$data5.ind.names)
grp2=score[score$PC1<=0,]
length(grp2$data5.ind.names)

#extraire les individus du data
g1data=data5[data5$ind.names %in% grp1$data5.ind.names]
g2data=data5[data5$ind.names %in% grp2$data5.ind.names]

#calculer frequence par groupe
g1freq=glMean(g1data)
g2freq=glMean(g2data)

g1f=data.frame(data5$loc.names, g1freq, 1-g1freq)
g2f=data.frame(data5$loc.names, g2freq, 1-g2freq)
names(g1f)=c("SNP", "q","p")
names(g2f)=c("SNP", "q","p")
delta_p=abs(g1f$p-g2f$p)

bestg1=g1f[g1f$SNP %in% best$data5.loc.names, ]
bestg2=g2f[g2f$SNP %in% best$data5.loc.names, ]
bestdelta_p=abs(bestg1$p-bestg2$p)

png("delta_p-golbal.png")
hist(delta_p, breaks=100)
dev.off()
png("delta_p-Q8.png")
hist(bestdelta_p, breaks=100)
dev.off()

# Refaire PCA sans les plus fort SNP

##newdata=data5[data5$loc.names %in% best$data5.loc.names] ##TODO marche pas avec data donc refaire vcf et output

perl -i.bak -pe "s/_[ACTG]//g" SNP-best0.02.txt 
vcftools --vcf pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf --exclude ../10-adegenet/SNP-best0.02.txt --out sansbest --recode
vcftools --vcf sansbest.recode.vcf --plink-tped --out sansbest
perl -i.bak -pe "s/^A\..-//g|s/-[0-9]+.A/\tA/g" sansbest.tfam 
plink --noweb --tfile sansbest --recode12 --out sansbest
plink --noweb --tfile sansbest --recodeA --out sansbest


#### PCA adegenet par pop

geno=read.table("final.gen", sep=" ", header=True)
pop=read.table("population.txt")

freq=char2genet(geno, pop$V1)
pcatot=dudi.pca(freq$tab, center = FALSE, scale = FALSE)
col=pcatot$li

write.table(pcatot$l1, "results-pca-pop.txt")
write.table(pcatot$c1, "results-pca-SNPpop.txt")

s.label(pcatot$li, xax = 6, yax = 7, sub = "5-6", lab = freq$pop.names)
dev.copy(pdf,"anguilla_acp_pop-axe5-6.pdf",width=7,height=5)
dev.off()
dev.off()

#PC1_2
graph_title="ACP in Anguilla anguilla"
n=nrow(col) # TO ANNOTATE THE SNP NUMBER ON GRAPH
x_title="PC1"
y_title="PC2"
graph_1<-ggplot(col,aes(x=col$Axis1 ,y=col$Axis2))
graph_1+geom_point()
graph_1+geom_point(aes(colour=freq$pop.names))+
scale_colour_manual(name="Populations",values=rainbow(16))+
labs(title=graph_title)+
labs(x=x_title)+
labs(y=y_title)

ggsave("anguilla_acp_PCA1_2.pdf",width=7,height=5,dpi=600)
dev.off()

### autre option pour analyse par pop, transformé fichier genepop avant l'analyse 
> data(microbov)
> obj <- genind2genpop(microbov,missing="chi2")
Converting data from a genind to a genpop object...
Replaced 0 missing values
...done.
> ca1 <- dudi.coa(as.data.frame(obj$tab),scannf=FALSE,nf=3)
> barplot(ca1$eig,main="Correspondance Analysis eigenvalues", col=heat.colors(length(ca1$eig)))

#### Find cluster adegenet

grp=find.clusters(data, max.n.clust=20)

### analyse descriminate

dapc1=dapc(data, eco$V6)

#### Fst bootsrap
fst=stamppFst(data, nboots = 500, percent = 95, nclusters = 3)
write.table(fst$Fsts, "fst-matrix-stampp.txt")
write.table(fst$Pvalues, "fst-pvalues-stampp.txt")
write.table(fst$Bootstraps, "fst-bootstrap-stampp.txt")

################ Arlequin #######################

../scripts/arlecore_linux/arlecore3513_64bit pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.arp amova.ars


##################### Bayescan ########################3
#Pour paire de pop
for i in *.geste; do ./BayeScan2.1_linux64bits $i -threads 16 -out_pilot -nbp 5 -pilot 1000 -burn 10000; done


### analyser output de bayescan avec script de bayescan
for i in *fst.txt; do cat $i > input_file.temp; R -q -e 'source("../../../scripts/BayeScan2.1/R functions/plot_R.r"); file="input_file.temp"; results=plot_bayescan(file,FDR=0.05); a=results$nb_outliers; write.table(a,"outliers.txt")'; mv Rplots.pdf $i.pdf; mv outliers.txt nb.outlier$i; done

touch all_comparisons_nb-outliers
for i in nb.outlier*; do echo -e "$i\t$(cut -d " " -f 2 $i)" >> all_comparisons_nb-outliers; done
rm nb.outlierpop*

### analyser output de bayescan avec script de thierry

for i in pop*fst.txt; do cat $i > input_file.temp; R -q -e 'source("../../../../08-script-bayescan.r")' ; mv recode-bayescan.txt recode-$i; mv Q6.txt Q6-$i; mv bayescan_anguilla-ecotype.png $i.png; done

touch all_comparisons_Q6.txt
for i in Q6*; do echo -e "$i\n$(grep -v X $i)" >> all_comparisons_Q6.txt; done
grep -v "Q6" all_comparisons_Q6.txt | sort| uniq -c| sort
perl -pe "s/\n/ /g|s/Q6/\nQ6/g" all_comparisons_Q6.txt | grep "scaffold4640-20255"

################ Plink #######################

### Association #######

# change 6eme colunm du tfam avec phenotype
plink --tfile pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17 --assoc --noweb --allow-no-sex

plink --tfile pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17 --fisher --noweb --allow-no-sex

plink --tfile pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17 --model --noweb --allow-no-sex

plink --tfile pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17 --logistic --noweb --allow-no-sex

################ Admixture #######################

for K in 1 2 3 4 5; do ../../../scripts/admixture --cv  $K | tee log${K}.out; done

#### Compare test d'association et acp

assoc=read.table("plink.assoc", header=T)
Q8=read.table("pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.Q8.txt", header=T)
delta=abs(assoc$F_A-assoc$F_U)
as=data.frame(assoc,delta)

bestassoc=as[as$P<0.01,]
write.table(bestassoc$SNP, "SNP-associes-0.005", row.names=F, quote=F)

aQ8=as[as$SNP %in% Q8$x,]
hist(aQ8$delta, proba=TRUE, col="gold",xlab="Delta p",ylab="Fréquence", breaks=20, main="Distribution du Delta p avec 450 SNP")

resacp=Q8$x[!(Q8$x %in% bestassoc$SNP)]
resas=bestassoc$SNP[!(bestassoc$SNP %in% Q8$x)]
hist(as$delta[as$SNP %in% resacp])
hist(as$delta[as$SNP %in% resas])

#### Recode vcf for phasing avec plink

head -n 9 pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf  > header.vcf 
grep -E -v "#" pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf | cut -f 3- > essai.recode.vcf
grep -E -v "#" pop1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-17.vcf | cut -f 1,2,3 | perl -pe "s/^scaffold//g|s/\|size[0-9]//g" | awk '{print $1+$2;}'| perl -pe "s/^/1\t/g"> col-position.txt 
paste col-position.txt essai.recode.vcf > 16poprecode.vcf
cat header.vcf 16poprecode.vcf > 16pop-recode.vcf
rm header.vcf
rm 16poprecode.vcf
rm col-position.txt
rm essai.recode.vcf

vcftools --vcf 16pop-recode.vcf --plink-tped --out 16pop
perl -i.bak -pe "s/^A\..-//g|s/-[0-9]+.A/\tA/g" 16pop.tfam
plink --noweb --tfile 16pop --recode12 --out 16pop
plink --noweb --tfile 16pop --recodeA --out 16pop

plink --tfile 16pop --hap-window 100 --noweb --allow-no-sex

