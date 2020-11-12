#!/bin/bash
# Create a population structure file for the ADMIXTURE supervised analysis by adding
# a dashed line to the known ancestry population strucutre file for every unique sample in the cohort

family_pedigree_file=$1

known_ancestry_pop_file=$2

admixture_pop_file=$3

cat $known_ancestry_pop_file > "${admixture_pop_file}"

number_of_samples=$(($(cut -f 1 $family_pedigree_file | wc -l) - 1526))

for i in $(eval echo {1..$number_of_samples})
do
  echo "$(echo '-'; cat $admixture_pop_file)" > $admixture_pop_file
done
