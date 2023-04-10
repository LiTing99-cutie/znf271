mkdir DAPARS2
mv Dapars2_* DAPARS2

# lncRNA cds annotation (run in Mac) later run in NanJing
cp /home/user/data2/lit/lncRNA_orf_prot_co_filter_genomic.txt annotation/

# znf271p hg38 gtf annotation (including CDS)
grep ENST00000399070.3 ../../02-APA/annotation/gencode.v41.annotation.gtf >ENST00000399070.3.before_anno.gtf
# edit in VScode editor according to /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/cds_lncRNA.txt
# ENST00000399070.3.reanno.gtf
cat ENST00000399070.3.before_anno.gtf ENST00000399070.3.reanno.gtf > ENST00000399070.3.gtf
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/ENST00000399070.3.gtf /home/user/data2/rbase/ucsc2/htdocs/data/lit

# ENSMMUG00000049532 proximal and distal annotation
grep ENSMMUT00000088941 /home/user/data/lit/database/in_house/rheMac10Plus/rheMac10Plus.addgeneName.gtf > ENSMMUT00000088941.gtf
cp ENSMMUT00000088941.gtf ENSMMUT00000088941.proximal.gtf
cp ENSMMUT00000088941.gtf ENSMMUT00000088941.distal.gtf
# cat ../../devo_compare/liftover/rheMac8.ZNF271.2Region.rheMac10Plus.gtf
# proximal start 47191612
# distal start 47189015

# ENSMMUG00000049532 proximal and distal annotation
zless /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf.gz | grep -w ENSMUST00000074941.7 > ENSMUST00000074941.7.gtf
cp ENSMUST00000074941.7.gtf ENSMUST00000074941.7.proximal.gtf
cp ENSMUST00000074941.7.gtf ENSMUST00000074941.7.distal.gtf
# less ../../devo_compare/mouse/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.gtf | grep -w Zfp35
# proximal end 24004729
# distal end 24005370

# rabbit gtf
mkdir rabbit && wget http://ftp.ensembl.org/pub/release-106/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.106.gtf.gz -P rabbit
zless /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rabbit/Oryctolagus_cuniculus.OryCun2.0.106.gtf.gz | grep -w ENSOCUT00000034050 > ENSOCUT00000034050.gtf
grep -w ENSOCUG00000029705 

# dog gtf
cd /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/
mkdir dog && wget http://ftp.ensembl.org/pub/release-106/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz  -P dog
cd /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/annotation
zless /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz | grep -w ENSCAFT00845012304 > ENSCAFT00845012304.gtf

# gtfToGenePred
gtfToGenePred dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.GenePred
gtfToGenePred -genePredExt dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.GenePredExt
genePredToBed dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.GenePred dog/Canis_lupus_familiaris.ROS_Cfam_1.0.wholeGene_annotation.bed

# IDmapping
cut -f 1,12 dog/Canis_lupus_familiaris.ROS_Cfam_1.0.106.GenePredExt > dog/IDmapping.txt
echo -e "#name\tname2" > dog/dapars2.idmapping.header.txt
cat dog/dapars2.idmapping.header.txt dog/IDmapping.txt > dog/IDmapping.wh.txt

# dog proximal and distal gtf
cp ENSCAFT00845012304.gtf ENSCAFT00845012304.proximal.gtf
cp ENSCAFT00845012304.gtf ENSCAFT00845012304.distal.gtf

# generate dapars annotation for human 
cd /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/
mkdir human_dapars_anno && cd human_dapars_anno
gtfToGenePred ../gencode.v41.basic.annotation.gtf gencode.v41.basic.annotation.GenePred
genePredToBed gencode.v41.basic.annotation.GenePred gencode.v41.basic.annotation.wholeGene_annotation.bed
cut -f 1,12 ../gencode.v41.basic.annotation.gpe > IDmapping.txt
cat ../dog/dapars2.idmapping.header.txt IDmapping.txt > IDmapping.wh.txt
