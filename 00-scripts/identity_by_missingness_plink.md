# Diagnose your missing genotypes pattern with PLINK's identity-by-missingness (IBM) analysis.

## Convert vcf file to plink tped
```
vcftools --vcf yourvcffile.vcf --plink-tped --out your_plinkfile
```

## Cluster on missing
```
plink --tfile your_plinkfile --cluster-missing --mds-plot 4 --out plink.ibm
```

## Add colums with needed info do .mds file (POP, lane, sexe...)
```
mds <- read.table ("plink_58pop.ibm.mds", header = TRUE)
pop_map <- read.table ("pop_map_all.tsv", header = TRUE)
mds_modified <- cbind (mds, pop_map$POP_ID)
strata.select = as.factor(pop_map$POP_ID) # factor you choose to look at the IBM, here the POP_ID
```

## Make graph in R
```
pIBM = ggplot(mds_modified, aes(x = C1, y = C2, colour = strata.select))
pIBM = p+geom_point() + labs(title = "MultiDimensional Scaling (MDS)\n of Identity by Missing (IBM) lake trout 58 pop") + 
  labs(x = "PC1") + labs(y = "PC2") + 
  theme(axis.title.x = element_text(size = 12,family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 12,family = "Helvetica", face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold")
        )
pIBM
ggsave("IBM_58_pop.pdf", width = 15, height = 15, dpi = 600, units = "cm", useDingbats = F)
```
