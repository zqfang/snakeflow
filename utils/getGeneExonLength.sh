cat GRCh38.gencode.v33.annotation.gtf | 
awk '{if($3=="exon") print $10"\t"$14"\t"$16"\t"$5-$4+1}' | 
sed -e 's/"//g' -e 's/;//g' | 
bedtools groupby -i - -g 1 -c 2 -o sum > Exon_lengths.txt

            