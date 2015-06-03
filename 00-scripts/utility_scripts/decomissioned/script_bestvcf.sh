#!/bin/bash

### Global variables

FILTERED_VCF=$1
BESTLOCI_FILE=$2
QVALUE=$3

### Best vcf

vcftools --vcf $FILTERED_VCF \
        --snps $BESTLOCI_FILE\
        --out $FILTERED_VCF-best$QVALUE \
        --recode


### Produce files for PLINK (tfam + tbed)

    vcftools --vcf $FILTERED_VCF-best$QVALUE.recode.vcf --plink-tped \
        --out $FILTERED_VCF-best$QVALUE

perl -i.bak -pe "s/^A\..-//g|s/-[0-9]+.A/\tA/g" $FILTERED_VCF-best$QVALUE.tfam 


### Modify tfam for structure populationnelle


### Produce file for adegenet
 	plink --noweb \
        --tfile $FILTERED_VCF-best$QVALUE \
        --recode12 \
        --out $FILTERED_VCF-best$QVALUE
	plink --noweb \
        --tfile $FILTERED_VCF-best$QVALUE \
        --recodeA \
        --out $FILTERED_VCF-best$QVALUE

### Produce files for Admixture
    # Remove loci
    plink --noweb \
        --tfile $FILTERED_VCF-best$QVALUE \
        --make-bed \
        --out $FILTERED_VCF-best$QVALUE
	
### Produce files PGD
	# Launch PGDSpider
	java -Xmx8g -Xms512M -jar ../../../scripts/PGDSpider2-cli.jar \
        -inputfile $FILTERED_VCF-best$QVALUE.recode.vcf\
        -inputformat VCF \
        -outputfile $FILTERED_VCF-best$QVALUE.pgd \
        -outputformat PGD \
        -spid ../../vcf_to_pgd.spid.ready

### Produce files for GenePop

    # For all individuals
    # Launch PGDSpider
    java -Xmx8g -Xms512M -jar ../../../scripts/PGDSpider2-cli.jar \
        -inputfile $FILTERED_VCF-best$QVALUE.pgd \
        -inputformat PGD \
        -outputfile $FILTERED_VCF-best$QVALUE.genepop \
        -outputformat GENEPOP \
        -spid ../../../00-input_files/pgd_to_genepop.spid


### Produce files for Arlequin

    # For all individuals
    # Launch PGDSpider
    java -Xmx8g -Xms512M -jar ../../../scripts/PGDSpider2-cli.jar \
        -inputfile $FILTERED_VCF-best$QVALUE.pgd \
        -inputformat PGD \
        -outputfile $FILTERED_VCF-best$QVALUE.arp \
        -outputformat ARLEQUIN \
        -spid ../../../00-input_files/pgd_to_arlequin.spid

rm *.vcfidx
