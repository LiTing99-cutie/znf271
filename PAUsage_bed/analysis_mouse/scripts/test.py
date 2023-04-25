import requests, sys
import pandas as pd
import os
from optparse import OptionParser

# Parse command line arguments
parser = OptionParser()
parser.add_option("-i", "--input", dest="input_file", help="gene_id to predict",default='/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict/to_predict_gene_id.txt')
parser.add_option("-o", "--output", dest="output_path", help="Output file path",default='/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict/cdna/')
parser.add_option("-a", "--archive", dest="archive_name", help="archive monthYear",default='Jul2022')
parser.add_option("-s", "--species", dest="species_name", help="species name or alias",default='homo_sapiens')
(options, args) = parser.parse_args()
print(options.archive_name)