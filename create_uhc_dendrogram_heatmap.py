from cospec import plot_uhc_heatmap
from cospec import plot_uhc_dendrogram
from figures import spectrum_map
import collections
from scipy.cluster.hierarchy import dendrogram
from fastcluster import linkage
from sklearn.metrics.pairwise import cosine_similarity

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac
import scipy.stats as sc
import glob
import csv
import argparse

#COSMIC_csv_reader.py
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#Author of cospec.py: linakim [at] mit.edu
#July 26th, 2018

SUBSTITUTION = 0
CONTEXT = 1
NORMALIZED_PROPORTION = 3

parser = argparse.ArgumentParser(description="""See the Juptyer Notebook for help""",
								 epilog="James was here")
parser.add_argument("-i", "--input",
					default="../../muts_csv/",
					help="The path to the directory containing select contexted sample .csv's. Make sure it ends with a '/'")
args = parser.parse_args()

cluster_names = []
spec_list = collections.OrderedDict()

for file in glob.glob(args.input+"*.csv"):
	filename = file.replace(args.input, "")
	cluster_names.append(filename.replace("_contexted", "").replace('.csv',''))
	file_dictionary = {}
	with open(file) as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if row[SUBSTITUTION] == "Substitution": continue
			subcon = "(" + str(row[SUBSTITUTION]) + " ," + str(row[CONTEXT]) + ")"
			file_dictionary[subcon] = float(row[NORMALIZED_PROPORTION])
	spec_list[file.replace("_contexted", "").replace(".csv","")] = file_dictionary

print("Plotting dendrogram...")
plot_uhc_dendrogram(spec_list, cluster_names)
print("Plotting heatmap. If this takes a while you probably are comparing too many samples.")
plot_uhc_heatmap(spec_list, cluster_names)
print("DONE")

