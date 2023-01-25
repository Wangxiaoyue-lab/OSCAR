# 1 fastqc
mkdir cleandata
mkdir fastqc
cat ./cleandata/files.txt | while read line
do
		#echo $line
		fastqc ./cleandata/$line --outdir ./fastqc -t 8
done



# 2 STAR
mkdir STAR
cat .cleandata/samples.txt | while read line
do
        echo $line
        mkdir ./STAR/$line
        STAR --genomeDir ./STAR/GRCm38_sjdbOverhang100/ --runThreadN 10 --readFilesIn ./cleandata/${line}_*_1.clean.fq.gz ./cleandata/${line}_*_2.clean.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix ./STAR/$line/$line_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard 
done



# 3 samtools
cat ./cleandata/samples.txt | while read line
do
        echo $line
        bampath=./STAR/$line/Aligned.sortedByCoord.out.bam
        #echo $bampath
        samtools sort $bampath -o ./STAR/$line/${line}_sort.bam -@ 20 -m 10G
        samtools index ./STAR/$line/${line}_sort.bam ./STAR/$line/${line}_sort.bam.bai
done



# 4 featureCounts
mkdir featurecounts
featureCounts -T 20 -a ./gencode.vM20.primary_assembly.annotation.gtf -o ./featurecounts/Featurecounts.txt -p ./STAR/*/*_sort.bam 2>./featurecounts/featurecounts.log