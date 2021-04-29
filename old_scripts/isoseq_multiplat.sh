#!/bin/bash

source $(dirname $0)/ath_env.conf
cd $outDir
mkdir -p $outDir/log

################### Functions 
checkFileSuffix(){
    filename=$(basename -- "$1")
    extension="${filename##*.}"
    if [ $extension == "$2" ];then
    	true
    else
    	false
    fi
}
################### End

inputDir=""
platform=$(echo "$platform"| tr a-z A-Z)

if [ $platform == "PACBIO" ];then	####### PacBio
	###### 1. Extract reads
	$pacPlat=$(echo "$pacPlat"| tr a-z A-Z)
	##### 1.1 Determine sequencing strategy
	if [ $pacPlat == "RSII" ];then
		# RS2
		find $inputDir -name '*.bax.h5' >input.fofn
		bax2bam -f input.fofn -o input	# Output input.subreads.bam
	elif [ $pacPlat == "SEQUEL" ];then
		# Sequel
		find $inputDir -name '*.subreads.bam' >subreadBam.fofn
		bamtools merge -list subreadBam.fofn -out input.subreads.bam
	else
		echo "You should input the right PacBio sequence strategy!!\n"
		exit 1
	fi
	##### 1.2 Extract raw reads
	echo -e "\n\n[INFO] $(date) Extract reads..."
	inputBam=input.subreads.bam
	
	samtools fastq --threads $thread $inputBam >input.rawsubreads.fastq 2>log/samtools.rawsubreads.fastq.log
	fqNameMapping.pl -p $sample/ input.rawsubreads.fastq >input.rawsubreads.fq 2>rawsubreads.idmapping
	echo -e "\n[INFO] $(date) Extract rawreads and subreads done"
	
	##### 1.3 Get CCS.fa
	echo -e "\n\n[INFO] $(date) Extract CCS..."
	ccs $inputBam CCS.bam -j $thread --noPolish --minPasses $minPass
	samtools view CCS.bam | cut -f1,13 | sed 's/np:i://' >name2pass.tsv
	samtools fastq --threads $thread CCS.bam 2>log/samtools.ccs.fastq.log | fqNameMapping.pl -p $sample/ >CCS.fq 2>CCS.idmapping
	awk 'NR%4==1{sub(/^@/,">",$0);print}NR%4==2{print}' CCS.fq >CCS.fa
	echo -e "\n[INFO] $(date) Extract CCS done"
	
	###### 2. Get full-length(FL) reads and pre-mapping evaluation
	##### 2.1 Get FL reads
	echo -e "\n\n[INFO] $(date) Get full-length (FL) reads..."
	hmmerWrapperCommonParams="--primer_search_window $primerSearchWin --min-score $primerSearchMinScore --cpus $thread --change-seqid --input_filename CCS.fa --output_filename FL.fa --directory hmmer/ --must-see-polyA"
	(time hmmer_wrapper.py $hmmerWrapperCommonParams --primer_filename $primer 2>log/hmmer_wrapper.log) 2>log/hmmer_wrapper.time
	summarize_primer_info.py FL.fa.primer_info.txt >log/summarize_primer_info.log
	mv FL.fa.primer_info.txt primer_info.tsv

	##### 2.2 Identify chimera
	(time chimera_finder.py --min_dist-from_end 50 --cpus $thread --primer_filename $primer --input_filename FL.fa --directory hmmer_chimera/ 2>log/chimera_finder.FL.log) 2>log/chimera_finder.time
	mv FL.fa.non_chimera.fa FLnoC.fa
	mv FL.fa.is_chimera.fa FLC.fa

	##### 2.3 Trim polyA
	faTrimPA.pl -w $windowSize -f $fraction FLnoC.fa 2>log/faTrimPA.log | faFilter -minSize=$minCcsLength stdin paTrimmed.fa
	(
		awk '$3<$2' log/faTrimPA.log | sed 's/^>//' | skyjoin - primer_info.tsv 2>/dev/null | awk '$10>$9{print $1"\t"$10-$9+$2-$3}'
		awk '$3<$2' log/faTrimPA.log | sed 's/^>//' | filter.pl -o /dev/stdin <(grep -v '^ID' primer_info.tsv) | awk '$7!="NA"&&$8!="NA"&&$8-$7>0{print $1"\t"$8-$7}'
	) | sort -n >paTailLength.tsv
	cut -f2 paTailLength.tsv | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null
	echo -e "\n[INFO] $(date) Get full-length (FL) reads done"
	
	#grep -v '^ID' primer_info.tsv | awk '$7!="NA"&&$8!="NA"&&$8-$7>0{print $1"\t"$8-$7}' | sort -n >paTailLength.tsv
	#cut -f2 paTailLength.tsv | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null
	#cp FLnoC.fa paTrimmed.fa
	##### 2.4 Pre-mapping evaluation
	echo -e "\n\n[INFO] $(date) Pre-mapping evaluation..."
	mkdir -p evaluation
	cd evaluation
	(
		faCount ../paTrimmed.fa >paTrimmed.faCount
		grep -v '^#' paTrimmed.faCount | grep -v total | awk '{print $1"\t"($4+$5)/$2*100}' >GC_of_reads.log
		cut -f2 GC_of_reads.log | distrCurve.R -d -m='GC Content of Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_reads.pdf 2>/dev/null

		GC_across_read.pl -i 20 ../paTrimmed.fa >GC_across_read.log
		cut -f2- GC_across_read.log | box.R -stack -nJ -ho=50 -m='GC Content across Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_read.pdf 2>/dev/null

		BQ_across_read.pl ../CCS.fq | tee BQ_across_read.log | cut -f2- | box.R -nJ -stack -m='Quality across Reads' -x=Interval -y=Quality -w=11 -p=BQ_across_read.pdf 2>/dev/null
	)
	echo -e "\n[INFO] $(date) Pre-mapping evaluation done"
	#lima --isoseq --dump-clips --no-pbi -j $thread $outPre.ccs.bam $primerFa $outPre.demux.bam
	#samtools fastq -0 -@ $thread $outPre.demux.5p--3p.bam > $outPre.fl.fq

	##### 3. Process with different correction tools
	if [ $anaTool == "proovread" ];then
		## 3.1 proovread pipeline
		echo -e "\n[INFO] $(date) Proovread correction..."
		mkdir -p proovread/{log,tmp}
		cd proovread
		SeqChunker -s $chunkSize -l 1 -o pb-%03d.fq ../paTrimmed.fa
		for i in `ls pb-*.fq`;do
			chunkName=$(echo "$i"| cut -d'.' -f1)
			proovread -l $i -s $leftReads -s $rightReads --prefix $chunkName -t $thread 2>log/$chunkName.proovread.log
			cd tmp/
			ln -sf ../$chunkName/*.trimmed.fa ./
			cd ../
		done
		rm pb-*.fq

		cat tmp/*.trimmed.fa > proovread.fa
		minimap2 -d ref.mm2 $refGenome
		minimap2 -ax splice:hq -uf ref.mm2 -t $thread proovread.fa > aln.sam
		#samtools view -bS aln.sam | samtools sort -@ $thread -o aln.sorted.bam
		#samtools index aln.sorted.bam
		echo -e "\n[INFO] $(date) Proovread correction done"
	elif [ $anaTool == "MECAT" ];then
		## 3.2 MECAT pipeline

		#ccs $inputBam $outPre.ccs.bam -j $thread --noPolish --minPasses $minPass
		#lima --isoseq --dump-clips --no-pbi --ccs -j $thread $outPre.ccs.bam $primerFa $outPre.demux.bam
		#isoseq3 refine $outPre.demux.F0--R0.bam $primerFa $outPre.flnc.bam

		#samtools fastq -0 -@ $thread $outPre.demux.F0--R0.bam > $outPre.fl.fq
		#samtools fastq -0 -@ $thread $outPre.flnc.bam > $outPre.flnc.mecat.fq

		#bamtools convert -format fasta -in $subreadsBam > m54148_170925_091709.subreads.fa
		#$pbtranscript classify --min_score 8 --summary summary.txt --report primer_info.csv --flnc FLnoC.fa --nfl notFL.fa --outDir hmmer --cpus 8 m54148_170925_091709.subreads.fa m54148_170925_091709.subreads.draft.fa --primer primer.fasta

		echo -e "\n[INFO] $(date) MECAT correction..."
		mkdir -p mecat
		cd mecat
		mecat2pw -j 0 -d ../paTrimmed.fa -o condidate.txt -w wrk_dir -t $thread -x $mecatMode -n $mecatChunkSize -a $minOverlap
		mecat2cns -i 0 -t $thread -x $mecatMode -p $partNum -r $minReadScore -a $minOverlap -c $minReadCov -l $minReadLength condidate.txt ../paTrimmed.fa corrected.fa

		minimap2 -d ref.mm2 $refGenome
		minimap2 -ax splice:hq -uf ref.mm2 corrected.fa -t $thread > aln.sam
		echo -e "\n[INFO] $(date) MECAT correction done"
		#samtools view -bS aln.sam | samtools sort -@ 8 -o aln.sorted.bam
		#samtools index aln.sorted.bam
	elif [ $anaTool == "isoseq3" ];then
		## 3.3 isoseq3 pipeline
		echo -e "\n[INFO] $(date) IsoSeq3 correction..."
		mkdir -p isoseq3/primer
		cd isoseq3
		cp $primerFa primer/
		makePrimerForLimaAndIsoseq.py $primerFa 1>primerIsoseq.fa 2>primerLima.fa
		lima --isoseq --dump-clips --no-pbi --ccs -j $thread CCS.bam primerLima.fa demux.bam 
		isoseq3 refine --require-polya $outPre.demux.5p--3p.bam primerIsoseq.fa flnc.bam 
		isoseq3 cluster -j $thread flnc.bam unpolish.flnc.bam
		isoseq3 polish -j $thread unpolish.flnc.bam $inputBam polished.flnc.bam

		samtools fastq -0 -@ $thread polished.flnc.bam > polished.flnc.fq
		minimap2 -d ref.mm2 $refGenome
		minimap2 -ax splice:hq -uf ref.mm2 -t $thread polished.flnc.fq > aln.sam
		echo -e "\n[INFO] $(date) IsoSeq3 correction done"
		#samtools view -bS aln.sam | samtools sort -@ $thread -o aln.sorted.bam
		#samtools index aln.sorted.bam
	fi
	samtools view -bS aln.sam | samtools sort -@ 8 -o aln.sorted.bam
	samtools index aln.sorted.bam
	cd ../
elif [ $platform == "NANOPORE" ];then	####### Nanopore
	inputFq=nanopore.fq
	if [ $anaTool == "proovread" ];then
		#ccs $inputBam $outPre.ccs.bam --noPolish --minPasses $minPass -j $thread
		#lima --isoseq --dump-clips --no-pbi -j $thread $outPre.ccs.bam $primerFa $outPre.demux.bam
		#isoseq3 refine $outPre.demux.F0--R0.bam $primerFa $outPre.flnc.bam

		#samtools fastq -0 -@ $thread $outPre.demux.F0--R0.bam > $outPre.fl.fq

		###### 4.1 proovread pipeline
		echo -e "\n[INFO] $(date) Proovread correction..."
		mkdir -p proovread/{log,tmp}
		cd proovread
		SeqChunker -s $chunkSize -l 1 -o pb-%03d.fq ../$inputFq
		for i in `ls pb-*.fq`;do
			chunkName=$(echo "$i"| cut -d'.' -f1)
			proovread -l $i -s $leftReads -s $rightReads --prefix $chunkName -t $thread 2>log/$chunkName.proovread.log
			cd tmp/
			ln -sf ../$chunkName/*.trimmed.fa ./
			cd ../
		done
		rm pb-*.fq

		cat tmp/*.trimmed.fa > proovread.fa
		ln -sf proovread.fa paTrimmed.fa
		
		minimap2 -d ref.mm2 $refGenome
		minimap2 -ax splice:hq -uf ref.mm2 -t $thread proovread.fq > aln.sam
		echo -e "\n[INFO] $(date) Proovread correction done"
		#samtools view -bS aln.sam | samtools sort -@ $thread -o aln.sorted.bam
		#samtools index aln.sorted.bam
	elif [ $anaTool == "MECAT" ];then
		#ccs $inputBam $outPre.ccs.bam --noPolish --minPasses $minPasses -j $thread
		#lima --isoseq --dump-clips --no-pbi --ccs -j $thread $outPre.ccs.bam $primerFa $outPre.demux.bam
		#isoseq3 refine $outPre.demux.F0--R0.bam $primerFa $outPre.flnc.bam

		#samtools fastq -0 -@ $thread $outPre.demux.F0--R0.bam > $outPre.fl.fq

		###### 4.2 mecat pipeline
		echo -e "\n[INFO] $(date) MECAT correction..."
		mkdir -p mecat
		cd mecat
		mecat2pw -j 0 -d ../$inputFq -o condidate.txt -w $outPre.wrk_dir -t $thread -x $mecatMode -n $mecatChunkSize -a $minOverlap 
		mecat2cns -i 0 -t $thread -x $mecatMode -p $partNum -r $minReadScore -a $minOverlap -c $minReadCov -l $minReadLength condidate.txt ../$inputFq corrected.fa
		ln -sf corrected.fa paTrimmed.fa
		
		minimap2 -d ref.mm2 $refGenome
		minimap2 -ax splice -uf -k14 ref.mm2 corrected.fa > aln.sam
		echo -e "\n[INFO] $(date) MECAT correction done"
		#samtools view -bS aln.sam | samtools sort -@ $thread -o aln.sorted.bam
		#samtools index aln.sorted.bam
	else
		echo "You should input the right tool name that you want to analysis the Nanopore sequence!!!"
		exit 1
	fi
	echo -e "\n\n[INFO] $(date) Nanopore post-correction evaluation..."
	mkdir -p evaluation
	cd evaluation
	(
		faCount paTrimmed.fa >paTrimmed.faCount
		grep -v '^#' paTrimmed.faCount | grep -v total | awk '{print $1"\t"($4+$5)/$2*100}' >GC_of_reads.log
		cut -f2 GC_of_reads.log | distrCurve.R -d -m='GC Content of Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_reads.pdf 2>/dev/null

		GC_across_read.pl -i 20 ../paTrimmed.fa >GC_across_read.log
		cut -f2- GC_across_read.log | box.R -stack -nJ -ho=50 -m='GC Content across Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_read.pdf 2>/dev/null

		BQ_across_read.pl ../CCS.fq | tee BQ_across_read.log | cut -f2- | box.R -nJ -stack -m='Quality across Reads' -x=Interval -y=Quality -w=11 -p=BQ_across_read.pdf 2>/dev/null
	)
	echo -e "\n[INFO] $(date) Nanopore post-correction evaluation done"
	
	samtools view -bS aln.sam | samtools sort -@ 8 -o aln.sorted.bam
	samtools index aln.sorted.bam
	cd ../
elif [ $platform == "NGS" ];then	####### RNA-seq
	rm -R RNA-seq 2>/dev/null
	mkdir -p RNA-seq/indexDir
	cd RNA-seq
	if checkFileSuffix $refAnno "gtf";then
		refAnnoGTF=$refAnno
	else
		refAnnoBase=$(basename -- "$refAnno")
		refAnnoName="${refAnnoBase%.*}"
		refAnnoGTF=${refAnnoName}.gtf
		gffread -T -o $refAnnoGTF $refAnno
	fi
	hisat2_extract_splice_sites.py $refAnnoGTF >${refAnnoName}.ss
	hisat2_extract_exons.py $refAnnoGTF >${refAnnoName}.exon
	hisat2_build -p $thread --ss ${refAnnoName}.ss --exon ${refAnnoName}.exon $refGenome indexDir/${refAnnoName}
	mappingName=$(echo $leftReads|cut -d'.' -f1)
	hisat2Para="-x indexDir/${refAnnoName} -1 $leftReads -2 $rightReads --dta -p $thread -S ${mappingName}.sam"
	hisat2 $hisat2Para 
	samtools sort -@ $thread -o ${mappingName}.sorted.bam ${mappingName}.sam
	samtools index ${mappingName}.sorted.bam
	cd ../
else 
	echo "You should input the right platform name that your data belongs to!!\n"
	exit 1
fi


##### Filter and refine alignments
mkdir -p mapping
cd mapping
echo -e "\n[INFO] $(date) Filter and refine alignments..."
### SpliceGrapher build classification module to identify correct splice site
build_classifiers.py -m $refAnnoGTF -f $refGenome -n 2000 -d gt,gc
samAddTag.pl --checkHardClip --coverage --identity --unmapped unmapped.sam aln.sam 2>lengthInconsistent.sam \
| samtools view -buS - | samtools sort -m 15G -@ $thread - >mapped.sorted.bam
sam2bed.pl -t CV,ID mapped.sorted.bam >mapped.bed12+

readsFilter.pl -r 0.8 mapped.bed12+ 2>discarded.bed12+ >uniq.bed12+
(
	samtools view -H mapped.sorted.bam
	samtools view mapped.sorted.bam | filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i
) | samtools view -buS - | sort -k 3,3 -k 4,4n - > uniq.sorted.sam
### SpliceGrapher use the classification module to filter sam file to get final correct sam file
sam_fileter.py uniq.sorted.sam classifier.zip -o filter.uniq.sorted.sam
samtools sort -m 8G -@ 5 filter.uniq.sorted.sam >uniq.sorted.bam

strandAdjust.pl -f $genomeSeq -g ../../raw.gpe -c 0.8 -s 2 uniq.bed12+ >strandAdjusted.bed12+
paste <(cut -f4,6 uniq.bed12+) <(cut -f6 strandAdjusted.bed12+) | awk '$2!=$3{print $1}' >strandConflict.readName
filter.pl -o strandConflict.readName -2 4 uniq.bed12+ >strandConfirm.bed12+

junctionScoringParams="-r ../../raw.gpe exonRealign.bed12+"
if [ "$junction" ];then
   removeMisAlignExon.pl -g ../../raw.gpe -j $junction strandConfirm.bed12+ | fillMissExon.pl -g ../../raw.gpe -b -j $junction >exonRealign.bed12+
   juncScoringParams="$juncScoringParams -j $junction"
else
   removeMisAlignExon.pl -g ../../raw.gpe -j /dev/null strandConfirm.bed12+ | fillMissExon.pl -g ../../raw.gpe -b -j /dev/null >exonRealign.bed12+
fi
juncConsensus.pl -s <(juncScoring.pl $juncScoringParams) exonRealign.bed12+ >consensus.bed12+

#cp consensus.bed12+ noGap.bed12+
bed2gpe.pl consensus.bed12+ | gpeFeature.pl --exon | filter.pl -o /dev/stdin -1 4 -2 4 consensus.bed12+ >noGap.bed12+

### This step can't remove internal priming and the identification of poly-A at 3' likely to be in-correct
dnaOrInternalPrimingContFilter.pl -b $twoBit -t <(sed 's/^@//' ../paTailLength.tsv) -r dnaOrInternalPrimingContFilter.log noGap.bed12+ >deCont.bed12+ 2>dnaCont.bed12+
(
	samtools view -F 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /(\d+)S$/ && $1 >30'
	samtools view -f 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /^(\d+)S/ && $1 >30'
) | filter.pl -o /dev/stdin -2 4 deCont.bed12+ >processed.bed12+
cut -f 4 processed.bed12+ | grep_fa_by_id.pl -i paTrimmed.fa -l - >processed.fa
cd ../
echo -e "\n[INFO] $(date) Filter and refine results done"


###### Remove artifact isoforms with machine learning
echo -e "\n\n[INFO] $(date) Remove artifact isoforms..."
##### The input is processed.fa and the output is artifact-removed fasta file
"Run sqanti pipeline"
mkdir sqanti
cd sqanti
bedToGenePred processed.bed12+ processed.gp
genePredToGtf file processed.gp processed.gtf
sqanti_qc.py processed.gtf $refAnnoGTF $refGenome -g -d ./ -o sqanti
sqanti_filter.py sqanti_classification.txt
grep_fa_by_id.pl -i processed.fa -l sqanti_classification.txt_curatedTranscriptome.txt >sqanti_artifact_removed.fa
cd ../
echo -e "\n[INFO] $(date) Remove artifact isoforms done"


###### Collapse the isoforms and map the corrected isoforms to the collapsed isoforms to quantify the isoforms. The input is the artifact-removed fasta file and output is the collapsed.fa and read count for each isoform(use the count to calculate the isoform expresse and PacBio can't be quantified)
mkdir collapse
cd collapse
ln -sf ../sqanti/sqanti_artifact_removed.fa
ln -sf ../sqanti/processed.bed12+
minimap2 -ax splice -uf -k14 ref.mm2 sqanti_artifact_removed.fa | sort -k 3,3 -k 4,4n - > artifact_removed.sorted.sam
collapse_out="tofu"

##collapse isoforms and generate $collapse_out.collapsed.gff, $collapse_out.collapsed.group.txt, $collpase_out.rep.fa which can be use in the following analysis
collapse_isoforms_by_sam.py --input sqanti_artifact_removed.fa -s artifact_removed.sorted.sam --dun-merge-5-shorter -o $collapse_out
echo -e "pbid\tcount_fl" > $collapse_out.collapsed.abundance.txt
paste <(awk '{print $1}' $collapse_out.collapsed.group.txt) <(awk -F , '{print NF}' $collapse_out.collapsed.group.txt ) >> $collapse_out.collapsed.abundance.txt
cut -f1 artifact_removed.sorted.sam > artifact_removed.read.lst
filter.pl -o $collapse_out.ignored_ids.txt artifact_removed.read.lst | seqtk subseq sqanti_artifact_removed.fa - > sqanti_artifact_removed.ignore_id_removed.fa
filter.pl -o $collapse_out.ignored_ids.txt processed.bed12+ > processed.ignore_id_removed.bed12+
cd ../
#minimap2 collpase.rep.fa sqanti_artifact_removed.fa  ### Count the number the read mapped to the sqanti_artifact_removed.fa and calculate expression

mdkir multi_sample_merge
# merge all identified read in each sample and generate new merged files: all_samples.chained.gff, all_samples.chained.rep.fq, all_samples.chained.rep.fa, all_samples.chained_ids.txt
chain_samples.py sample.config norm_nfl --dun-merge-5-shorter

################################# Above steps correct the isoforms and get the hq_isoform.fa to carry out following analysis


###### 5. Post-mapping Evaluation
echo -e "\n\n[INFO] $(date) Post-mapping evaluation..."
mkdir -p evaluation/length
cd evaluation
##### 5.2 Coverage across Genes with Differential Length
ln -sf ../mapping/processed.ignore_id_removed.bed12+
ln -sf ../collapse/sqanti_artifact_removed.ignore_id_removed.fa
faCount sqanti_artifact_removed.ignore_id_removed.fa > sqanti_artifact_removed.ignore_id_removed.faCount
grep -v '^#' sqanti_artifact_removed.ignore_id_removed.faCount | grep -v total | cut -f2 >length/Read.lst
gpe2bed.pl ../../raw.gpe | bedLength.pl | cut -f13 >length/ReferenceGene.lst
cut -f1-12 processed.ignore_id_removed.bed12+ | bedLength.pl | cut -f13 >length/processed.ignore_id_removed.lst
distrCurves.R -x1=0 -x2=10000 -d -x='Binned Length (limited in 0-10000)' -w=15 length/*.lst -b=150 -p=LengthDistribution.curve.pdf 2>/dev/null
boxes.R -ng -no length/*.lst -p=LengthDistribution.box.pdf 2>/dev/null

geneCoverage.pl -g ../../raw.gpe processed.ignore_id_removed.bed12+ >coverage.log
awk '$3<20000' coverage.log | select.pl --index 3,2 | binBox.R -b=10 -w=15 -m='Coverage across Genes' -x='RefSeq Gene Length (<20000)' -y=Coverage -p=coverage.pdf 2>/dev/null
##### 5.3 Distance to TSS/TTS across Reads with Different Length
distanceToEnd.pl -g ../../raw.gpe -f sqanti_artifact_removed.ignore_id_removed.faCount processed.ignore_id_removed.bed12+ | awk '$2!="NA"' >distance2Gene.log
select.pl --index 4,2 distance2Gene.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TSS' -x='Binned Read Length' -y=Distance -p=distance2TSS.pdf 2>/dev/null
select.pl --index 4,3 distance2Gene.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TTS' -x='Binned Read Length' -y=Distance -p=distance2TTS.pdf 2>/dev/null
cut -f2 distance2Gene.log >DistanceToTSS.dis
cut -f3 distance2Gene.log >DistanceToTTS.dis
distrCurves.R *.dis -x1=-100 -x2=100 -b=10 -w=15 -d -x=Distance -y='Density of Reads' -m='Distance to TSS/TTS' -p=distance2End.pdf 2>/dev/null
distrCurve.R <DistanceToTTS.dis -d -x='Distance (limited in -500 to 500)' -y='Density of Reads' -m='Distance to TTS' -x1=-500 -x2=500 -b=10 -p=distance2EndOfTTS.pdf 2>/dev/null
rm *.dis
##### 5.4 Expression in different genome region
readsAssigner.pl -g ../../raw.gpe -s processed.ignore_id_removed.bed12+ >assign.bed12+
function genomeAbundance(){
	exonicReadN=$(awk '$15=="E"' $1 | wc -l)
	intronicReadN=$(awk '$15=="I"' $1 | wc -l)
	intergenicReadN=$(awk '$15=="IG"' $1 | wc -l)
	totalReadN=$(wc -l <$1)
	exomeLength=$(gpeFeature.pl -e $2 | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}')
	intromeLength=$(gpeFeature.pl -i $2 | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}')
	genomeLength=$(awk 'BEGIN{SUM=0}{SUM+=$2}END{print SUM}' $3)
	genomeLength=$(($genomeLength * 2))	### why multiple 2
	intergenomeLength=$(($genomeLength - $exomeLength - $intromeLength))
	exomeRPKM=$(echo "scale=3;$exonicReadN / ($exomeLength/1000) / ($totalReadN / 10^6)" | bc)
	intromeRPKM=$(echo "scale=3;$intronicReadN / ($intromeLength/1000) / ($totalReadN / 10^6)" | bc)
	intergenomeRPKM=$(echo "scale=3;$intergenicReadN / ($intergenomeLength/1000) / ($totalReadN / 10^6)" | bc)
	echo -e "Exome\t$exomeRPKM"
	echo -e "Introme\t$intromeRPKM"
	echo -e "Intergenome\t$intergenomeRPKM"
}
genomeAbundance assign.bed12+ ../../raw.gpe $chrSize >genomeAbundance.log
bar.R <genomeAbundance.log -anno -c=darkgreen -f=white -m=genomeAbundance -x='Genomic Context' -y=RPKM -p=genomeAbundance.pdf 2>/dev/null
##### 5.5 Saturation of Sequencing
totalReadN=$(wc -l <processed.ignore_id_removed.bed12+)
function saturation(){
	readsAssigner.pl -g $1 $2 >assign.multLine.bed12+
	transN=$(wc -l <$1)
	geneN=$(cut -f12 $1 | sort -u | wc -l)
	step1=$(($totalReadN/80))
	step2=$(($totalReadN/40))
	halfOfTotal=$(($totalReadN/2))
	printf '' >trans.discover
	printf '' >gene.discover
	printf '' >exon.discover
	for(( i=0; i<=$halfOfTotal; i=$(($i+$step1)) ));do
		shuf $2 | head -n $i >shuf.bed12+
		transReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $16}' | sort -u | wc -l)
		geneReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $17}' | sort -u | wc -l)
		discoverRate=$(echo "scale=3;$transReads/$transN" | bc)
		echo -e "$i\t$discoverRate" >>trans.discover
		discoverRate=$(echo "scale=3;$geneReads/$geneN" | bc)
		echo -e "$i\t$discoverRate" >>gene.discover
		printf "$i\t" >>exon.discover
		exonDiscoverRate.pl -g $1 shuf.bed12+ >>exon.discover
	done
	for((; i<=$totalReadN; i=$(($i+$step2)) ));do
		shuf $2 | head -n $i >shuf.bed12+
		transReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $16}' | sort -u | wc -l)
		geneReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $17}' | sort -u | wc -l)
		discoverRate=$(echo "scale=3;$transReads/$transN" | bc)
		echo -e "$i\t$discoverRate" >>trans.discover
		discoverRate=$(echo "scale=3;$geneReads/$geneN" | bc)
		echo -e "$i\t$discoverRate" >>gene.discover
		printf "$i\t" >>exon.discover
		exonDiscoverRate.pl -g $1 shuf.bed12+ >>exon.discover
	done
	rm -f shuf.bed12+
}
saturation ../../raw.gpe processed.ignore_id_removed.bed12+
lines.R *.discover -m='Saturation of Sequencing' -x='Reads Count' -y='Discovery Rate' -w=15 -lgPos=top -p=discoveryRate.pdf 2>/dev/null
##### 5.6 Sequencing Error Rate
awk 'BEGIN{OFS="\t"}{print $4,100-$13*$14/100}' processed.ignore_id_removed.bed12+ | tee errorRate.log | cut -f2 | distrCurve.R -d -m='Rate of Mismatch & Indel for Mapping' -x=Rate -y='Fraction of Reads' -p=errorRate.pdf 2>/dev/null
##### 5.7 Read Count per Gene
geneRPKM.pl -g ../../raw.gpe processed.ignore_id_removed.bed12+ | tee RPKM.bed6+ | cut -f7 >readCount.lst		######The method calculating the RPKM may have some mistakes
distrCurve.R <readCount.lst -m='Gene Count at Read Count' -x='Read Count' -y='Gene Count' -xl=10 -yl=10 -b=500 -ng -p=readCountPerGene.pdf 2>/dev/null
cd ../
echo -e "\n[INFO] $(date) Post-mapping evaluation done"


###### 6. PA Analysis
echo -e "\n\n[INFO] $(date) PA Analysis..."
##### 6.1 PA identification
awk 'BEGIN{FS=OFS="\t"}$6=="+"{$2=$3-1}$6=="-"{$3=$2+1}{print $1,$2,$3,$6}' processed.ignore_id_removed.bed12+ | sort -u | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"ClevageSite"NR,0,$4}' >cleavage.bed6
bedClosest.pl cleavage.bed6 | tee adjacentCleavage.tsv | cut -f13 | distrCurve.R -xl=2 -d -x='Binned Distance between Adjacent Cleavage Site' -p=adjacentCleavageDistr.pdf 2>/dev/null		###This line can be remove for the adjacent is not important!!
rm -f paCluster.count
for i in $(seq 5 1 100) $(seq 100 10 10000) $(seq 10000 1000 100000);do
	echo -e $i"\t"$(paCluster.pl -d $i -m mode -w 3 processed.ignore_id_removed.bed12+ | wc -l) >>paCluster.count
done
line.R <paCluster.count -x=Distance -y=Count -p=paClusterCount.pdf 2>/dev/null
paCluster.pl -d $paDistance -m mode -w 3 processed.ignore_id_removed.bed12+ | tee paCluster.bed8+ | awk '$5>1{print $3-$2}' | distrCurve.R -d -x='Cluster Size (limited in 1-100)' -y='Density' -m='Distribution of Cluster Size' -x1=0 -x2=100 -b=1 -p=paClusterSize.pdf
awk 'BEGIN{FS=OFS="\t"}{$2=$7;$3=$8;print}' paCluster.bed8+ | cut -f1-6,9- >PA.bed6+

PAClassbyRead.pl -a ../evaluation/assign.bed12+ <(cut -f1-8 paCluster.bed8+) >paCluster.type.bed8+ 2>singleExonReadWithExonInMEread.bed12+		##"PAClassbyRead.pl" line 89-92 may be incorrect
awk '$9!="SE"' paCluster.type.bed8+ >paCluster.deSingleExonRead.bed8+
awk '$9=="SE"' paCluster.type.bed8+ >paCluster.singleExonRead.bed8+
##### 6.2 Motif around PA
mkdir -p motif
up=100
down=100
for seq in A T C G;do
	awk '$5>1' PA.bed6+ | motifAroundPA.pl -u $up -d $down -m $seq -b $twoBit >motif/$seq.nucleotide 2>/dev/null
done
lines.R -w=15 -y=Frequency -x='Distance Relative to PA' -m='Distribution of Nucleotide Frequency around PA' -p=motif/nucleotide.pdf motif/*.nucleotide 2>/dev/null

up2=50
down2=0
for seq in AATAAA AAATAA ATAAAA ATTAAA ATAAAT TAATAA ATAAAG AAAATA CAATAA ATAAAC AAAAAA AAAAAG;do
	awk '$5>1' PA.bed6+ | motifAroundPA.pl -u $up2 -d $down2 -m $seq -b $twoBit >motif/$seq.PAS 2>/dev/null
done
lines.R -p=motif/PAS.1.pdf motif/{AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null
lines.R -p=motif/PAS.2.pdf motif/{ATAAAG,AAAATA,CAATAA,ATAAAC,AAAAAA,AAAAAG}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null
paste motif/ATAAAC.PAS <(cut -f2 motif/ATAAAG.PAS) <(cut -f2 motif/CAATAA.PAS) <(cut -f2 motif/AAAAAA.PAS) <(cut -f2 motif/AAAAAG.PAS) <(cut -f2 motif/AAAATA.PAS) \
| awk '{print $1"\t"$2+$3+$4+$5+$6+$7}' >motif/Other.PAS
lines.R -p=motif/PAS.pdf motif/{Other,AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null
echo -e "\n[INFO] $(date) PA Analysis done"





###### 8. Identify ASEs and DE transcripts
##### 8.1 Identify alternative splicing events
echo -e "\n\n[INFO] $(date) Alternative splicing events identifying..."
mkdir -p ASE
cd ASE
readGroup.pl -g ../raw.gpe ../evaluation/assign.bed12+ >unambi.bed12+ 2>ambiguous.bed12+
cut -f1-12,21 unambi.bed12+ | tee readGrouped.bed12+ | filter.pl -o <(cut -f4 ../PA/singleExonReadWithExonInMEread.bed12+; cut -f4 ../PA/paCluster.singleExonRead.bed8+ | tr ',' "\n") -2 4 >deSingleExonRead.bed12+		###Remove single exon reads
3endRevise.pl -p ../PA/PA.bed6+ deSingleExonRead.bed12+ | tee deSingleExonRead.3endRevised.bed12+ | paGroup.pl >paGrouped.tsv 2>paGrouped.bed6
mkdir -p PB
getSingleIr.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/IR.bed6+
geneWithIR=$(cut -f7 PB/IR.bed6+ | sort -u | wc -l)
deSingleExonReadCount=$(cut -f13 deSingleExonRead.bed12+ | sort -u | wc -l)
geneWithIrRatio=$(awk 'BEGIN{printf "%.4f", '$geneWithIR/$deSingleExonReadCount'}')
getSE.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/SE.bed12+
getA5SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A5SS.bed6+
getA3SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A3SS.bed6+
awk '$9>1 && $11>1 && $5>=100 && $5<=900{print $3-$2}' PB/A5SS.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
awk '$3-$2>=20 && $9>1 && $11>1{print $5/1000}' PB/A5SS.bed6+ | hist.R -x='Incusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
ln -sf SE.bed12+ PB/SE.confident.bed12+
ln -sf A5SS.bed6+ PB/A5SS.confident.bed6+
ln -sf A3SS.bed6+ PB/A3SS.confident.bed6+
awk '{print $3-$2}' PB/A5SS.confident.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
awk '$9>1 && $11>1{print $5/1000}' PB/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
awk '{print $3-$2}' PB/A3SS.confident.bed6+ | hist.R -p=PB/A3SS.size.pdf 2>/dev/null
awk '$9>1 && $11>1{print $5/1000}' PB/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A3SS.InclusionRatio.pdf 2>/dev/null
echo -e "\n[INFO] $(date) Alternative splicing events identifying done"

##### 8.2 Alternative splicing events characterization
echo -e "\n\n[INFO] $(date) APE Characterization..."
mkdir -p APE/characterization
cd APE/characterization

#### Get event
awk '{print $1":"$4":"$6}' ../PB/SE.bed12+              | grep -v Novel >PB.SE.lst
awk '{print $1":"$4":"$6}' ../PB/A5SS.bed6+             | grep -v Novel >PB.A5SS.lst
awk '{print $1":"$4":"$6}' ../PB/A3SS.bed6+             | grep -v Novel >PB.A3SS.lst
awk '{print $1":"$4":"$6}' ../PB/IR.bed6+               | grep -v Novel >PB.IR.lst
awk 'BEGIN{FS="\t";OFS=":"}{split($4,array,",")}length(array)>1&&$3!~/^chr/{print $1,$2,$3,$4}' ../paGrouped.tsv >PB.APA.lst
venn.R PB.IR.lst -p=IR.pdf 2>/dev/null
venn.R PB.APA.lst -p=APA.pdf 2>/dev/null
ln -sf PB.SE.lst confident.SE.lst
ln -sf PB.A5SS.lst confident.A5SS.lst
ln -sf PB.A3SS.lst confident.A3SS.lst

#### Classify APE into annotated or novel according to reference gene model
awk '{print $0"\t"$12}' ../../../raw.gpe | gpe2bed.pl -p >reference.bed12+
getSingleIr.pl -g ../../../raw.gpe reference.bed12+ | tee IR.reference.bed6+ | awk '{print $1":"$4":"$6}' >IR.reference.lst
paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6
paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv
getSE.pl -g ../../../raw.gpe reference.bed12+ | tee SE.reference.bed12+ | awk '{print $1":"$4":"$6}' >SE.reference.lst
getA5SS.pl -g ../../../raw.gpe reference.bed12+ | tee A5SS.reference.bed6+ | awk '{print $1":"$4":"$6}' >A5SS.reference.lst
getA3SS.pl -g ../../../raw.gpe reference.bed12+ | tee A3SS.reference.bed6+ | awk '{print $1":"$4":"$6}' >A3SS.reference.lst

filter.pl -o IR.reference.lst PB.IR.lst >IR.novel.lst
filter.pl -o IR.reference.lst PB.IR.lst -m i >IR.anno.lst
#awk '{split($4,array,",")}length(array)>1' paGroup.novel.tsv
#awk '{split($4,array,",")}length(array)>1' paGroup.anno.tsv
filter.pl -o SE.reference.lst   confident.SE.lst        >SE.novel.lst
filter.pl -o SE.reference.lst   confident.SE.lst   -m i >SE.anno.lst
filter.pl -o A5SS.reference.lst confident.A5SS.lst      >A5SS.novel.lst
filter.pl -o A5SS.reference.lst confident.A5SS.lst -m i >A5SS.anno.lst
filter.pl -o A3SS.reference.lst confident.A3SS.lst      >A3SS.novel.lst
filter.pl -o A3SS.reference.lst confident.A3SS.lst -m i >A3SS.anno.lst

#### Get Statistic of APE
(echo -e "##AS ID is composed of Gene:Retained intron start-Retained intron end"
echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
awk '{print $1":"$4":"$6"\t"$0}' ../PB/IR.bed6+ | filter.pl -o IR.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
awk '{print $1":"$4":"$6"\t"$0}' ../PB/IR.bed6+ | filter.pl -o IR.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
) >statistic.IR.tsv
(echo -e "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
awk 'BEGIN{OFS="\t"}{split($4,array,",")}length(array)>1{print $1,$2,"Known",$3,$4}' paGroup.anno.tsv
awk 'BEGIN{OFS="\t"}{split($4,array,",")}length(array)>1{print $1,$2,"Novel",$3,$4}' paGroup.novel.tsv
) >statistic.APA.tsv
(echo -e "##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start"
echo -e "##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]"
echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene"
awk '{print $1":"$4":"$6"\t"$0}' ../PB/SE.confident.bed12+ | filter.pl -o SE.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$14}'
awk '{print $1":"$4":"$6"\t"$0}' ../PB/SE.confident.bed12+ | filter.pl -o SE.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$14}'
) >statistic.SE.tsv
(echo -e "##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end"
echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
awk '{print $1":"$4":"$6"\t"$0}' ../PB/A5SS.confident.bed6+ | filter.pl -o A5SS.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
awk '{print $1":"$4":"$6"\t"$0}' ../PB/A5SS.confident.bed6+ | filter.pl -o A5SS.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
) >statistic.A5SS.tsv
(echo -e "##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end"
echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
awk '{print $1":"$4":"$6"\t"$0}' ../PB/A3SS.confident.bed6+ | filter.pl -o A3SS.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
awk '{print $1":"$4":"$6"\t"$0}' ../PB/A3SS.confident.bed6+ | filter.pl -o A3SS.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
) >statistic.A3SS.tsv

#### Event Decompose and Motif Evaluation
awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' PB.IR.lst    >IR.splicesite
awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' IR.anno.lst  >IR.anno.splicesite
awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' IR.novel.lst >IR.novel.splicesite
seDecompose.pl        confident.SE.lst   >SE.inc.splicesite   2>SE.exc.splicesite
anssDecompose.pl -n 5 confident.A5SS.lst >A5SS.inc.splicesite 2>A5SS.exc.splicesite
anssDecompose.pl -n 3 confident.A3SS.lst >A3SS.inc.splicesite 2>A3SS.exc.splicesite
splicesite2figure

#### Distance of PA to TTS
awk '$4!~/^chr/' ../paGrouped.bed6 | awk 'BEGIN{OFS="\t"}$6=="+"{print $1,$3-1,$3,$4,$5,$6}$6=="-"{print $1,$2,$2+1,$4,$5,$6}' >pbPA.bed6
bedtools closest -a <(sort -k1,1 -k2,2n pbPA.bed6) -b <(gpeFeature.pl --tts ../../../raw.gpe|sort -k1,1 -k2,2n) -s -D a | select.pl -i 13,4 | sort -u | tee pbPA2TTS.tsv | cut -f1 | box.R -ng -nJ -no -y='Distance to TTS' -p=pbPA2TTS.pdf
echo -e "\n[INFO] $(date) APE Characterization done"
echo -e "\n\nAll analysis done"

##### 8.3 Identify differential alternative spliced/expressed isoforms between samples
strategy="TGS"
if [ $strategy == "TGS" ];then
	if [ $platform == "PACBIO" ];then
		echo "It can't be quantified"
	elif [ $platform == "NANOPORE" ];then
		echo "It can be quantified"
	fi
elif [ $strategy == "HYBRID" ];then
	if [ $platform == "PACBIO" ];then
		echo "It's PACBIO hybrid strategy"
	elif [ $platform == "NANOPORE" ];then
		echo "It's NANOPORE hybrid strategy"
	fi
elif [ $strategy == "NGS" ];then
	echo "It's NGS strategy"
fi


###### 9. Isoform visualization

