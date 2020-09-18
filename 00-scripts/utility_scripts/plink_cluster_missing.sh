#!/bin/bash
plink --vcf input_renamed.vcf \
    --cluster missing \
    --out input_renamed \
    --mds-plot 4 \
    --allow-extra-chr
