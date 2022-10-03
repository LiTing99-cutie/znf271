#!/bin/sh
usage(){
  cat << EOF
Description: Plot heatmap using deeptools plotHeatmap.
Author: Li Yumei, 2018/3/6
Options:
  -i  FILE    The input matrix to plot heatmap with header and row names.
  -u  INT     Upstream length of the plot
  -d  INT     Downstream length of the plot
  -b  INT     Body length of the plot
  -o  FILE    The output file pdf name.
  -a  STRING  The parameters passed to plotHeatmap except -m and -out.
  -h          Show this help information
EOF
    exit 0
}

while getopts "hi:u:d:b:o:a:" OPTION
do
    case $OPTION in
        h) usage;;
        i) input=$OPTARG;;
		u) upstream=$OPTARG;;
		d) downstream=$OPTARG;;
		b) body=$OPTARG;;
        o) outputPdf=$OPTARG;;
        a) heatmapArgs=$OPTARG;;
        ?) usage;;
    esac
done
shift $((OPTIND - 1))
colN=$(awk -v OFS="\t" '{print NF}' $input|head -n1|cut -f1)
rowN=$(wc -l $input|cut -f1 -d ' ')
echo "@{\"upstream\":$upstream,\"downstream\":$downstream,\"body\":$body,\"bin size\":1,\"verbose\":true,\"bin avg type\":\"mean\",\"missing data as zero\":false,\"min threshold\":null,\"max threshold\":null,\"scale\":1,\"skip zeros\":false,\"nan after end\":false,\"proc number\":20,\"unscaled 5 prime\":0,\"unscaled 3 prime\":0,\"sample_labels\":[\"\"],\"group_labels\":[\"\"],\"sample_boundaries\":[0,$colN],\"group_boundaries\":[0,$rowN]}" >tmp.header
awk -v OFS="\t" '{print "chr1","0","1",".",".",".",$0}' $input >tmp.input
cat tmp.header tmp.input >tmp.matrix
gzip tmp.matrix
plotHeatmap -m tmp.matrix.gz -out $outputPdf $heatmapArgs
rm tmp.matrix.gz tmp.header tmp.input
