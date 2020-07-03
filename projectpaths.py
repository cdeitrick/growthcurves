
from pathlib import Path
from typing import *

import pandas
from loguru import logger


import utilities


EXPECTED_FORMAT = "[strain].[consition].[plate].[replicate]"
AUC_COLUMN = 'auc_e'



class Filenames:
	""" Holds the filenames for important files/figures/etc. Also implements the methods to save tables/figures.

		Output File Stucture

	"""

	def __init__(self, folder: Path):
		self.table_format = '.tsv'
		self.figure_format = '.pdf'
		# TODO: Make sure these options are actually
		folder = utilities.checkdir(folder)
		self.folder_data = utilities.checkdir(folder / "data")
		self.folder_figure = utilities.checkdir(folder / "figures")

		self.strains = "WT,A244T,N274Y,N455K,P421L,tRNA".split(',')
		self.treatments = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')

		# Tables
		"""
			The statistics table describes the logistic population model each sample was fitted to.
			Columns:
			- `condition`: str
			- `strain`: str
			- `condition:strain`: str (can probably be removed since the `condition` and `strain` columns already exist.
			- `plate`: str
			- `replicate`: int
			- `N`: float
			- `k`: float
			- `r`: float
			- `sigma`: float
			- `auc_l`: float (The area under the logistic curve model.
			- `auc_e`: float (The area under the observed population from the plate reader.)
		"""

		# Basically the same as the `populationmodel` table, but with the sample id as well. Should remove the above table since it is redundant.
		self.filename_table_population = self.folder_data / ("populationmodel" + self.table_format)
		# Summarizes the regresson after applying an ordinary least-squares model to the input data.
		self.filename_table_regression_model = self.folder_data / "regression.txt"
		# Summarizes the results from the ANOVA analysis.
		self.filename_table_anova = self.folder_data / ("anova" + self.table_format)

		# Contains all paired tukey calulations. Tukey operates as a pairwise calculation of the difference in means for each variable pair.
		self.folder_tables_tukey = utilities.checkdir(self.folder_data / "tukey")
		self.filename_table_tukey = self.folder_tables_tukey / ("tukey" + self.table_format)
		self.filename_data_tukey = self.folder_data / "tukeyhsdresults.json"
		# Add a table with the maximum observed growth for each sample, aorted in ascending order. This should help
		# identify samples with little to no growth
		self.filename_table_maximum_growth = self.folder_data / "maximumgrowth.txt"
		self.filename_table_auc_statistics = self.folder_data / "auc_statistics.tsv"
		self.filename_table_growthcurve_models = self.folder_data / "growthcurve.model.tsv"

		# Figures
		self.folder_figures_growthcurves = utilities.checkdir(self.folder_figure / "growthcurves")
		self.folder_figures_tukey = utilities.checkdir(self.folder_figure / "tukey")
		self.folder_figures_anova = utilities.checkdir(self.folder_figure / "anova")
		self.filename_figure_sigmas = self.folder_figure / ("sigmas" + self.figure_format)
		self.filename_figure_qq = self.folder_figure / ("qq" + self.figure_format)
		self.filename_figure_anova_plot_main = self.folder_figure / ("anovaplot.main" + self.figure_format)
		self.filename_figure_anova_plot_groups = self.folder_figure / ("anovaplot.panel" + self.figure_format)

		self.filename_table_info = self.folder_data / "tableinfo.txt"
		self.controls = {
			'condition': 'RKS',
			'strain':    'WT'
		}



