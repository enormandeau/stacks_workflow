#!/bin/bash
# TODO remove unneeded files as needed

### Global variables
FILTERED_SUMSTATS=$(basename $(echo $1))
FILTERED_VCF=$(basename $(echo $2))
POPULATION_FILE=$(basename $(echo $3))
POPS_GROUPING=$4
POPS=$(echo $POPS_GROUPING | perl -pe 's/\-[0-9]+//g')
GROUPING=$(echo $POPS_GROUPING | perl -pe 's/[0-9]+\-//g')

# Input and output folders
INPUT_FOLDER=00-input_files
OUTPUT_LOCI=01-loci
OUTPUT_GROUPS=02-groups
OUTPUT_VCF=04-vcf
OUTPUT_FST=05-fst
OUTPUT_PLINK=06-admixture_plink
OUTPUT_GENEPOP=07-genepop
OUTPUT_STRUCTURE=08-structure
OUTPUT_BAYESCAN=09-bayescan
OUTPUT_ADEGENET=10-adegenet

# Create output folders
mkdir $OUTPUT_LOCI $OUTPUT_GROUPS $OUTPUT_VCF \
    $OUTPUT_BAYESCAN $OUTPUT_FST $OUTPUT_PLINK \
    $OUTPUT_GENEPOP $OUTPUT_STRUCTURE \
    $OUTPUT_ADEGENET 2> /dev/null


# Create group files
rm ./$OUTPUT_GROUPS/* 2> /dev/null

for g in $GROUPING
do
    touch ./$OUTPUT_GROUPS/$FILTERED_SUMSTATS-group_$g.txt
done

# Put individuals from populations into the appropriate group files
for g in $POPS_GROUPING
do
    pop=$(echo $g | cut -d "-" -f 1)
    group=$(echo $g | cut -d "-" -f 2)
    grep -E "[[:space:]]$pop$" ./$INPUT_FOLDER/$POPULATION_FILE | \
        cut -f 1 | \
        perl -sape 's/$/\t$g/' -- -g=$group >> \
        ./$OUTPUT_GROUPS/$FILTERED_SUMSTATS-group_$group.txt
done

cat $OUTPUT_GROUPS/$FILTERED_SUMSTATS-group* > \
    $OUTPUT_GROUPS/$FILTERED_SUMSTATS.txt

cut -f 1 $OUTPUT_GROUPS/$FILTERED_SUMSTATS.txt > \
    $OUTPUT_GROUPS/$FILTERED_SUMSTATS-1_column.txt

perl -pe "s/^A\..-//g|s/-[0-9]+\t/\t/g" $OUTPUT_GROUPS/$FILTERED_SUMSTATS.txt | \
    sort | \
    uniq > \
    $OUTPUT_GROUPS/$FILTERED_SUMSTATS.fst.pop.id.txt

# Replace '.' by './.'
if [ -f $INPUT_FOLDER/$FILTERED_VCF.diploid ]
then
    echo "diploid file exists"
else
    perl -pe "s/[[:space:]]\.:/\t.\/.:/g" $INPUT_FOLDER/$FILTERED_VCF > \
    $INPUT_FOLDER/$FILTERED_VCF.diploid
fi

### Filter vcf with best loci
cut -f 3,4 $INPUT_FOLDER/$FILTERED_SUMSTATS| uniq > \
    $OUTPUT_LOCI/$FILTERED_SUMSTATS-filtered-locis.txt 
cut -f 2 $INPUT_FOLDER/$FILTERED_SUMSTATS | \
    uniq > $OUTPUT_LOCI/$FILTERED_SUMSTATS-filtered-locis-2.txt 

vcftools --vcf $INPUT_FOLDER/$FILTERED_VCF.diploid \
    --positions $OUTPUT_LOCI/$FILTERED_SUMSTATS-filtered-locis.txt \
    --snps $OUTPUT_LOCI/$FILTERED_SUMSTATS-filtered-locis-2.txt \
	--keep $OUTPUT_GROUPS/$FILTERED_SUMSTATS-1_column.txt \
    --out $OUTPUT_VCF/$FILTERED_SUMSTATS \
    --recode

rm $OUTPUT_VCF/*.log 2> /dev/null

### Recode vcf file 

head -n 9 $OUTPUT_VCF/$FILTERED_SUMSTATS.recode.vcf > $OUTPUT_VCF/header.vcf 
grep -E -v "#" $OUTPUT_VCF/$FILTERED_SUMSTATS.recode.vcf > $OUTPUT_VCF/recode.vcf
awk '{print $1,$2}' $OUTPUT_VCF/recode.vcf| perl -pe "s/\|.+ | /-/g" > $OUTPUT_VCF/new-col
paste $OUTPUT_VCF/new-col $OUTPUT_VCF/recode.vcf > $OUTPUT_VCF/recode1.vcf
awk 'BEGIN {OFS="\t"} {$4=$1;print}' $OUTPUT_VCF/recode1.vcf > $OUTPUT_VCF/recode2.vcf
cut -f 2- $OUTPUT_VCF/recode2.vcf > $OUTPUT_VCF/final.vcf
cat $OUTPUT_VCF/header.vcf $OUTPUT_VCF/final.vcf > $OUTPUT_VCF/$FILTERED_SUMSTATS.vcf
grep -v -E "##"  $OUTPUT_VCF/$FILTERED_SUMSTATS.vcf | cut -f 3 > \
    $OUTPUT_LOCI/$FILTERED_SUMSTATS-keep-locis.txt

rm $OUTPUT_VCF/header.vcf
rm $OUTPUT_VCF/new-col
rm $OUTPUT_VCF/recode.vcf
rm $OUTPUT_VCF/recode1.vcf
rm $OUTPUT_VCF/recode2.vcf
rm $OUTPUT_VCF/final.vcf

### Produce files for PLINK (tfam + tbed)
vcftools --vcf $OUTPUT_VCF/$FILTERED_SUMSTATS.vcf --plink-tped \
    --out $OUTPUT_PLINK/$FILTERED_SUMSTATS

# TODO Modify .tfam file, second column must contain only population name
perl -i.bak -pe "s/^A\..-//g|s/-[0-9]+.A/\tA/g" $OUTPUT_PLINK/$FILTERED_SUMSTATS.tfam 

### Modify tfam for structure populationnelle
### Produce file for adegenete
plink --noweb \
    --tfile $OUTPUT_PLINK/$FILTERED_SUMSTATS \
    --recode12 \
    --out $OUTPUT_ADEGENET/$FILTERED_SUMSTATS
plink --noweb \
    --tfile $OUTPUT_PLINK/$FILTERED_SUMSTATS \
    --recodeA \
    --out $OUTPUT_ADEGENET/$FILTERED_SUMSTATS

### Produce files for Admixture
    # Remove loci
plink --noweb \
    --tfile $OUTPUT_PLINK/$FILTERED_SUMSTATS \
    --make-bed \
    --out $OUTPUT_PLINK/$FILTERED_SUMSTATS

### Produce files PGD
filename=$OUTPUT_GROUPS/$FILTERED_SUMSTATS.txt
perl -sape 's/INPUT_FILE/$f/' -- -f=$filename $INPUT_FOLDER/vcf_to_pgd.spid \
    > $INPUT_FOLDER/vcf_to_pgd.spid.ready

# Launch PGDSpider
java -Xmx8g -Xms512M -jar scripts/PGDSpider2-cli.jar \
    -inputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.vcf \
    -inputformat VCF \
    -outputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.pgd \
    -outputformat PGD \
    -spid $INPUT_FOLDER/vcf_to_pgd.spid.ready

### Produce files for Arlequin
# For all individuals
# Launch PGDSpider
java -Xmx8g -Xms512M -jar scripts/PGDSpider2-cli.jar \
    -inputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.pgd \
    -inputformat PGD \
    -outputfile $OUTPUT_GENEPOP/$FILTERED_SUMSTATS.genepop \
    -outputformat GENEPOP \
    -spid $INPUT_FOLDER/pgd_to_genepop.spid

### Produce files for GenePop
# For all individuals
# Launch PGDSpider
java -Xmx8g -Xms512M -jar scripts/PGDSpider2-cli.jar \
    -inputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.pgd \
    -inputformat PGD \
    -outputfile $OUTPUT_GENEPOP/$FILTERED_SUMSTATS.arp \
    -outputformat ARLEQUIN \
    -spid $INPUT_FOLDER/pgd_to_arlequin.spid
	

### Produce files for Structure
# For all individuals
# Launch PGDSpider
java -Xmx8g -Xms512M -jar scripts/PGDSpider2-cli.jar \
    -inputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.pgd \
    -inputformat PGD \
    -outputfile $OUTPUT_STRUCTURE/$FILTERED_SUMSTATS.structure \
    -outputformat STRUCTURE \
    -spid $INPUT_FOLDER/pgd_to_structure.spid

 cut -f 1,2 $OUTPUT_STRUCTURE/$FILTERED_SUMSTATS.structure | \
     grep -E "^\w" | \
     uniq > $OUTPUT_GROUPS/$FILTERED_SUMSTATS-genepop-id.txt 

### Produce files for Bayescan
# For all individuals
# Launch PGDSpider
java -Xmx8g -Xms512M -jar scripts/PGDSpider2-cli.jar \
    -inputfile $OUTPUT_VCF/$FILTERED_SUMSTATS.pgd\
    -inputformat PGD \
    -outputfile $OUTPUT_BAYESCAN/$FILTERED_SUMSTATS.geste\
    -outputformat GESTE \
    -spid $INPUT_FOLDER/pgd_to_bayescan.spid

head -n 4 $OUTPUT_BAYESCAN/$FILTERED_SUMSTATS.geste | \
    perl -pe "s/ons\]=16/ons\]=2/g" > header.txt 
perl -pe "s/\[pop\]\=.+/pop/g|s/^/c/g|s/\n//g|s/pop/\npop/g" $OUTPUT_BAYESCAN/$FILTERED_SUMSTATS.geste | \
    grep -E -v "\[|\]"  > geste.temp
possible_groups=$(cat $OUTPUT_GROUPS/$FILTERED_SUMSTATS-genepop-id.txt | cut -f 2 | uniq)

for g1 in $possible_groups; do for g2 in $possible_groups
do
    if [[ $g1 < $g2 ]]
    then
        awk -v i=$g1 -v j=$g2 'NR==i || NR==j' geste.temp | \
        perl -pe "s/c/\n/g|s/\n+\n/\na/g" | \
        perl -pe "s/apop/\n\[pop\]\=2/g|s/^pop/\[pop\]\=1/g|s/a//g" > \
        pop$g1-$g2-comp
    fi
done
done

for i in pop*-comp; do cat header.txt $i > $OUTPUT_BAYESCAN/$i.geste; done 
rm pop*-comp
rm geste.temp
rm header.txt

# Cleanup
rm $OUTPUT_VCF/*.vcfidx 2> /dev/null
rm $INPUT_FOLDER/*.vcfidx

