#!/bin/bash
set -e

# define usage
display_usage() { 
  echo -e "\nUsage: count_blocks.sh <merged.vcf.gz> <regions.txt> \n\nResults are printed to stdout\n"
} 

# if less than two arguments supplied, display usage 
if [  $# -le 1 ] 
then 
  display_usage
  exit 1
fi 
 
# check whether user had supplied -h or --help . If yes display usage 
if [[ ( $# == "--help") ||  $# == "-h" ]] 
then 
  display_usage
  exit 0
fi 

# loop thru regions list
while read region
do
  # get WT sibling hets
  WThet=`bcftools query -r $region \
  -i 'GT[0]="0/1"' \
  -f '%CHROM \t %POS \t %REF \t %ALT \t GT:[ %GT] \n' \
  $1 | wc -l | cut -f1 -d" "`
  
  # get Mut hets
  Muthet=`bcftools query -r $region \
  -i 'GT[1]="0/1"' \
  -f '%CHROM \t %POS \t %REF \t %ALT \t GT:[ %GT] \n' \
  $1 | wc -l | cut -f1 -d" "`
  
  # get Mut-specific homozygous mutations
  Muthomo=`bcftools query -r $region \
  -i 'GT[0]!="1/1" && GT[1]="1/1"' \
  -f '%CHROM \t %POS \t %REF \t %ALT \t GT:[ %GT] \n' \
  $1 | wc -l | cut -f1 -d" "`

  echo $region$'\t'$WThet$'\t'$Muthet$'\t'$Muthomo
done < $2
