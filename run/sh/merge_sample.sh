#!/usr/bin/sh

################################################
#File Name: ../merge_sample.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 30 Sep 2022 01:54:33 PM CST
################################################

set -eou pipefail

pushd /home/user/data2/lit/project/ZNF271/02-APA-1/analysis

# generate pattern file for grep 

for stage in $(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt|sed 's/ /_/g'|cut -f7|tail -n +2|sort|uniq)
do
less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt | \
sed 's/ /_/g' | awk -v stage=$stage '$7==stage' | awk '{print $2}' |sort| uniq > pattern/$stage.txt
ls target_bam/*znf271.bam|egrep -f pattern/$stage.txt|egrep -i "brain|cerebellum" > bam_lst/$stage.bam.lst 
done 

for stage in $(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt|sed 's/ /_/g'|cut -f9|tail -n +2|sort|uniq)
do
less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt | \
sed 's/ /_/g' | awk -v stage=$stage '$9==stage' | awk '{print $2}' |sort| uniq > pattern/$stage.txt
ls target_bam/*znf271.bam|egrep -f pattern/$stage.txt|egrep -i "brain|cerebellum" > bam_lst/$stage.bam.lst 
done 


for bam_lst in $(ls bam_lst/*bam.lst);do
stage=$(basename $bam_lst | sed 's/.bam.lst//g')
samtools merge -f -o bam_lst/$stage.forebrain.bam -b $bam_lst 
done

# bam to bedgraph
	for bam_lst in $(ls bam_lst/*bam.lst);do
	stage=$(basename $bam_lst | sed 's/.bam.lst//g')
	genomeCoverageBed -bg -trackline -ibam bam_lst/$stage.forebrain.bam -split > bedgraph/$stage.forebrain.bedgragh &
	done
	wait
#


# bedgraph normalization
	for bam_lst in $(ls bam_lst/*bam.lst);do
	stage=$(basename $bam_lst | sed 's/.bam.lst//g')
	m_r=$(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt| sed 's/ /_/g'|\
	grep $stage|cut -f 2)
	less bedgraph/$stage.forebrain.bedgragh|tail -n +2|\
	awk -v m_r=$m_r '{print $1,$2,$3,$4/m_r}' > bedgraph/$stage.forebrain.nm.bedgragh
	done
#


# bedgraph to bigwig
	for bam_lst in $(ls bam_lst/*bam.lst);do
	stage=$(basename $bam_lst | sed 's/.bam.lst//g')
	bedGraphToBigWig bedgraph/$stage.forebrain.nm.bedgragh /home/user/data/lit/database/public/genome/hg38/hg38_comChr.chrom.sizes \
	bedgraph/$stage.forebrain.nm.bw
	done
#



# 2022/12/23 adust bin_size

	# get bedgraph (bin size of 50 bp)
		for bam_lst in $(ls bam_lst/*bam.lst);do
		stage=$(basename $bam_lst | sed 's/.bam.lst//g')
		samtools index bam_lst/$stage.forebrain.bam && bamCoverage -b bam_lst/$stage.forebrain.bam -o bedgraph/$stage.forebrain.dt.bedgragh -of bedgraph -bs 50 -p 10 &
		done
		wait
	#

	# normalization (rpkm with bin size of 50 bp)
		for bam_lst in $(ls bam_lst/*bam.lst);do
		stage=$(basename $bam_lst | sed 's/.bam.lst//g')
		m_r=$(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt| sed 's/ /_/g'|\
		grep $stage|cut -f 2)
		less bedgraph/$stage.forebrain.dt.bedgragh|egrep -v 'alt|random|chrUn'|\
		awk -v m_r=$m_r '{print $1,$2,$3,$4/m_r/5*10^2}' > bedgraph/$stage.forebrain.dt.nm.bedgragh
		done

		for bam_lst in $(ls bam_lst/*bam.lst);do
		stage=$(basename $bam_lst | sed 's/.bam.lst//g')
		m_r=$(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt| sed 's/ /_/g'|\
		grep $stage|cut -f 2)
		less bedgraph/$stage.forebrain.dt.bedgragh|egrep -v 'alt|random|chrUn'|\
		awk -v m_r=$m_r '{print $1,$2,$3,$4/m_r}' > bedgraph/$stage.forebrain.dt.nm.rpm.bedgragh
		done

		# before birth sample:51
		# 0-2.57506 rpm
		# 0-51.5012 rpkm bin_size=50bp

		# after birth sample:28
		# 0-1.7823 rpm
		# 0-35.6461 rpkm bin_size=50bp

		# ajust the library size of merged samples in run/R/metadata.Rs

		# test bin size 
			for bin_size in 50 60 100 120 150 200 300;do
				bamCoverage -b bam_lst/After_birth.forebrain.bam -o bedgraph/After_birth.forebrain.dt.bedgragh -of bedgraph -bs $bin_size -p 10 2>/dev/null
				less bedgraph/After_birth.forebrain.dt.bedgragh | awk -v bs=$bin_size '{print $1,$2,$3,$4/1081/bs*10^3}'|cut -f 4 -d ' '|sort -k1,1n |tail -n1
			done

			# 56.9288
			# 46.4847
			# 34.4866
			# 28.4921
			# 26.9319
			# 25.7539
			# 16.5772
	#

	# get bigwig and cp to myhub
		for bam_lst in $(ls bam_lst/*bam.lst);do
		stage=$(basename $bam_lst | sed 's/.bam.lst//g')
		bedGraphToBigWig bedgraph/$stage.forebrain.dt.nm.bedgragh /home/user/data/lit/database/public/genome/hg38/hg38_comChr.chrom.sizes \
		bedgraph/$stage.forebrain.dt.nm.bw && \
		cp /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/bedgraph/$stage.forebrain.dt.nm.bw /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/rna_seq/develop_brain/znf271
		done
	
# 2022/12/23 run
	for bam_lst in $(ls bam_lst/*bam.lst);do
	stage=$(basename $bam_lst | sed 's/.bam.lst//g')
	m_r=$(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt| sed 's/ /_/g'|grep $stage|cut -f 2)
	bamCoverage -b bam_lst/$stage.forebrain.bam -o bedgraph/$stage.forebrain.dt.bedgragh -of bedgraph -bs 100 -p 10 2>/dev/null&& \
	less bedgraph/$stage.forebrain.dt.bedgragh|egrep -v 'alt|random|chrUn'|awk -v m_r=$m_r '{print $1,$2,$3,$4/m_r/100*10^3}' > bedgraph/$stage.forebrain.dt.nm.bedgragh && \
	bedGraphToBigWig bedgraph/$stage.forebrain.dt.nm.bedgragh /home/user/data/lit/database/public/genome/hg38/hg38_comChr.chrom.sizes bedgraph/$stage.forebrain.dt.nm.bw && \
	cp -f /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/bedgraph/$stage.forebrain.dt.nm.bw /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/rna_seq/develop_brain/znf271 &
	done