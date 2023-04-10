ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*bam | \
grep -f <(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt | awk '$2>0.885' | cut -f1 | tail -n +2) > frag_fil_develop_bam.lst

