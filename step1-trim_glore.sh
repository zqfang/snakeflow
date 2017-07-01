#!bin/bash
CleanPath="./clean_fastq"
CleanFastqcPATH="./qc/fastqc_clean"
R1suffix="R1.fastq.gz"
R2suffix="R2.fastq.gz"

mkdir -p $CleanFastqcPATH
mkdir -p logs/trim-galore
mkdir -p $CleanPath
ls raw_fastq/*R1.fastq.gz | while read var;
do
	id=${var%%_R*}
	sample=${id##*/}
        trim_galore -q 20 --phred33 --fastqc --stringency 1 --fastqc_args "--outdir $CleanFastqcPath" -e 0.1  --length 50 -o $CleanPath --paired ${id}_${R1suffix} ${id}_${R2suffix} > logs/trim-galore/${sample}.trimgalore.log
done


