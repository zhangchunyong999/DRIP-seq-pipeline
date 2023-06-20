dir=/beegfs/home/shilei/DRIPseq

mkdir {rawdata,qc,align,rmdupbam,genebody,bw}

cd $dir/rawdata

ls *gz|cut -d"_" -f 1|sort -u|while read id;do (trim_galore  -q 25 --phred33 --stringency 3 --length 30 -e 0.1 --paired ${id}_1.fq.gz ${id}_2.fq.gz --gzip -o $dir/cleandata/ ); done

cd $dir/qc

ls ../cleandata/*gz|xargs fastqc -t 20 -o ./

multiqc ./

cd $dir/cleandata/

bowtie2_index="/beegfs/home/shilei/ChIP-Seq/reference/hg19/hg19 "

ls *gz|cut -d"_" -f 1|sort -u|while read id;do  ls -lh ${id}_1_val_1.fq.gz ${id}_2_val_2.fq.gz;bowtie2 -p 20 -x $bowtie2_index -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz --no-unal|samtools sort  -O bam  -@ 40 -o - > ${id}.bam;done

mv *.bam $dir/align/

cd $dir/align/

ls  *.bam  | while read id ;do (nohup samtools flagstat $id > $(basename $id ".bam").stat & );done

ls *.bam|while read id; do ( picard MarkDuplicates ASO=coordinate REMOVE_DUPLICATES=true I=$(basename $id '.bam').bam O=$(basename $id '.bam').rmdup.bam M=$(basename $id '.bam').rmdup.txt ); done

mv *.rmdup.bam $dir/rmdupbam/

mv *.rmdup.txt $dir/rmdupbam/

cd $dir/rmdupbam/

ls  *.rmdup.bam  |xargs -i samtools index {} 

ls  *.rmdup.bam  | while read id ;do (nohup samtools flagstat $id > $(basename $id ".bam").stat & );done

ls *.rmdup.bam |while read id;do ( bamCoverage -p 20 --normalizeUsing RPKM -bs 10 -b $id -o ${id%%.*}.bw ); done

mv *.bw $dir/bw/

cd $dir/genebody

computeMatrix scale-regions  -p 20  \
--beforeRegionStartLength 25000 \
--regionBodyLength 25000 \
--afterRegionStartLength 25000 \
-R /beegfs/home/shilei/ChIP-Seq/reference/UCSC_Main_on_Human__knownGene.bed \
-S ../bw/*bw  \
--missingDataAsZero \
--skipZeros  -o matrix_body_25k.gz  \
--outFileSortedRegions regions_body_25k.bed

plotProfile -m matrix_body_25k.gz --perGroup -out sumprofile25kall.pdf --plotFileFormat pdf


