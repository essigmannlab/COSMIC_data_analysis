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

SUBSTITUTION = 0
CONTEXT = 1
NORMALIZED_PROPORTION = 3

cluster_names = []
spec_list = collections.OrderedDict()

for file in glob.glob("../../muts_csv/*.csv"):
	filename = file.replace("../../muts_csv/", "")
	cluster_names.append(filename.replace("_contexted", "").replace('.csv',''))
	file_dictionary = {}
	with open(file) as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if row[SUBSTITUTION] == "Substitution": continue
			subcon = "(" + str(row[SUBSTITUTION]) + " ," + str(row[CONTEXT]) + ")"
			file_dictionary[subcon] = float(row[NORMALIZED_PROPORTION])
	spec_list[file.replace("_contexted", "").replace(".csv","")] = file_dictionary

plot_uhc_dendrogram(spec_list, cluster_names)
plot_uhc_heatmap(spec_list, cluster_names)

