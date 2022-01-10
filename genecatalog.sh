for f in *.fa; do anvi-gen-contigs-database -T 4 -f $f -n ${f} -o ${f}_out.db; done
for f in *.db; do anvi-run-ncbi-cogs -c $f -T 4 --sensitive; done
for f in *.db; do anvi-run-kegg-kofams -c $f -T 4; done
for f in *.db; do anvi-run-pfams -c $f -T 4; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_COG20Category.txt --annotation-sources COG20_CATEGORY; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_COG20Function.txt --annotation-sources COG20_FUNCTION; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_KOfam.txt --annotation-sources KOfam; done
for f in *.db; do anvi-export-functions -c $f -o ${f}_Pfam.txt --annotation-sources Pfam; done
for f in *.db; do anvi-export-gene-calls -c $f --gene-caller prodigal -o ${f}_AllGeneCalls.txt; done

for f in *.db; do anvi-run-hmms -c $f -T 4; done
for f in *.db; do anvi-run-scg-taxonomy -c $f -T 4; done

for f in *KOfam.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_KOfamEvalue_filtered.txt;done
for f in *Pfam.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_PfamEvalue_filtered.txt;done
for f in *COG20Function.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_COG20FunctionEvalue_filtered.txt;done
for f in *COG20Category.txt; do awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' $f >${f}_COG20CategoryEvalue_filtered.txt;done


for f in *KOfamEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}KOfamCounts.txt; done
for f in *PfamEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}PfamCounts.txt; done
for f in *COG20FunctionEvalue_filtered.txt do cut -f4 $f| sort | uniq -c > ${f}COG20Function.txt; done
for f in *COG20CategoryEvalue_filtered.txt; do cut -f4 $f| sort | uniq -c > ${f}COG20Category.txt; done

for f in *KOfamCounts.txt; do sed 's/^ *//g' $f > ${f}KOfamCountsReformatted.txt;done


for f in *KOfamCountsReformatted.txt; do sed 's/ /;/' $f > ${f}KOfamCountsReformattedTwo.txt; done
 
for f in *KOfamCountsReformattedTwo.txt; do awk -F ";" '{n=split($1,a,";");for (i=1;i<=n;i++) print $2"\t"a[i]}' $f > ${f}FinalCounts.tsv; done

set -e

uniq_genes_file="uniq_genes.txt"

output_file="combined-output.tsv"


cat <( printf "function\n" ) ${uniq_genes_file} > ${output_file}


IFS=$'\n'


for curr_file in $(ls *FinalCounts.tsv); do


curr_sample_name=$(echo ${curr_file} | sed 's/.tsv//')


printf "${curr_sample_name}\n" > building-col.tmp


for curr_gene in $(cat ${uniq_genes_file}); do


if grep -m 1 -F -q "${curr_gene}" ${curr_file}; then

grep -m 1 -F "${curr_gene}" ${curr_file} | cut -f 2 >> building-col.tmp

else
 
  printf "0\n" >> building-col.tmp

if

done


paste ${output_file} building-col.tmp > building-output.tmp
mv building-output.tmp ${output_file}
rm building-col.tmp

done 

find . -name "*KOfamCounts.txt" -type f -delete
find . -name "*KOfamEvalue_filtered.txt" -type f -delete
find . -name "*Reformatted.txt" -type f -delete
find . -name "*ReformattedTwo.txt" -type f -delete
