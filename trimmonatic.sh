set -e
mkdir -p fastq_trim
ls fastq/*R1.fastq.gz | while read var;
    do
	id=${var%%_R*}
	id2=${id##*/}
	#echo -e $id2
	trimmomatic PE -threads 12 -phred33  ${id}_R1.fastq.gz ${id}_R2.fastq.gz  fastq_trim/${id2}_R1_paried.trim.fastq.gz fastq_trim/${id2}_R1_unpaired.trim.fastq.gz fastq_trim/${id2}_R2_paired.trim.fastq.gz fastq_trim/${id2}_R2_unpaired.trim.fastq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:50 | tee fastq_trim/${id2}.trim.log

	done
