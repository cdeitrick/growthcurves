import json
from io import StringIO
from pathlib import Path
from typing import *

import pandas
from loguru import logger
from statsmodels.regression import linear_model

import utilities
from graphics import AnovaPanelPlot, AnovaPlotNested, PlotGrowthcurves, other
from projectpaths import Filenames


class CleanTukey:
	def __init__(self, strains: Union[str, List[str]] = None, treatments: Union[str, List[str]] = None):
		"""
			Cleans up the tukey table.
		Parameters
		----------
		strains: Union[str,List[str]]
			Can be either a list or a comma-delimited string of strains.
		treatments: Union[str,List[str]]
			Can be either a list or a comma-delimited string of treatments.
		"""
		self.strains = "WT,A244T,N274Y,N455K,P421L,tRNA".split(',')
		self.treatments = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')

	# self.strains = strains.split(',') if isinstance(strains, str) else strains
	# self.treatments = treatments.split(',') if isinstance(treatments, str) else treatments

	@staticmethod
	def get_possible_combinations_single(iterable: List[Any]) -> List[Tuple[Any, Any]]:
		return [(i, j) for i in iterable for j in iterable if i != j]

	def add_missing_keys(self, table: pandas.DataFrame, name: str, labels: List[str]):
		label_1 = 'group1'
		label_2 = 'group2'

		current_table = table[table['name'] == name]
		current_labels = [(i, j) for i, j in zip(current_table[label_1].tolist(), current_table[label_2].tolist())]

		possible_labels = self.get_possible_combinations_single(labels)
		missing_labels = [i for i in possible_labels if i not in current_labels]
		debug = False
		if debug:
			logger.debug(f"Current labels: {current_labels[:5]}")
			logger.debug(f"Possible labels: {possible_labels[:5]}")
			logger.debug(f"Missing labels: {missing_labels[:5]}")
			logger.debug(f"There are {len(possible_labels)} possible combinations for '{name}'")
			logger.debug(f"There are {len(current_labels)} combinations in the table.")
			logger.debug(f"There are {len(missing_labels)} missing labels in the table.")

		reversed_table = current_table.copy()
		reversed_table[label_1] = current_table[label_2]
		reversed_table[label_2] = current_table[label_1]

		df = pandas.concat([current_table, reversed_table])

		return df

	def _prepare_tukey_results(self, tukey_table: pandas.DataFrame) -> pandas.DataFrame:
		""" Adds a few extra columns to the full tukey table to make future lookups easier."""

		def get_strain(s):
			if '-' in s:
				r = s.split('-')[1]
			elif s in self.strains:
				r = s
			else:
				r = None
			return r

		def get_treatment(s):
			if '-' in s:
				r = s.split('-')[0]
			elif s in self.treatments:
				r = s
			else:
				r = None
			return r

		reversed_table = tukey_table.copy()
		reversed_table['group1'] = tukey_table['group2']
		reversed_table['group2'] = tukey_table['group1']
		fulltable = pandas.concat([tukey_table, reversed_table])

		for column in ['group1', 'group2']:
			fulltable[column] = fulltable[column].apply(lambda s: s.strip())

		fulltable['strain1'] = fulltable['group1'].apply(get_strain)
		fulltable['strain2'] = fulltable['group2'].apply(get_strain)
		fulltable['condition1'] = fulltable['group1'].apply(get_treatment)
		fulltable['condition2'] = fulltable['group2'].apply(get_treatment)
		# Correct the 'reject' column. It's currently saved as a text string with additional whitespace.
		fulltable['reject'] = fulltable['reject'].apply(lambda s: s.strip() == 'True')
		return fulltable

	def clean(self, tukey_table: pandas.DataFrame) -> pandas.DataFrame:
		string_columns = ['group1', 'group2', 'name']
		for column in string_columns:
			tukey_table[column] = tukey_table[column].apply(lambda s: s.strip())

		group_plates = [f"plate{i}" for i in range(1, 13)]
		group_strains = self.strains
		group_conditions = self.treatments
		group_combined = [f"{i}-{j}" for i in group_strains for j in group_conditions]

		table_plate = self.add_missing_keys(tukey_table, 'plate', group_plates)
		table_strain = self.add_missing_keys(tukey_table, 'strain', group_plates)
		table_condition = self.add_missing_keys(tukey_table, 'condition', group_plates)
		table_combined = self.add_missing_keys(tukey_table, 'condition_strain', group_combined)

		fulltable = pandas.concat([table_plate, table_strain, table_condition, table_combined])

		return fulltable


class FigureWorkflow:
	""" Generates figures using the data from `GrowthCurveAnalysis."""

	def __init__(self, folder: Path, label_order: List[str] = None, groups: List[str] = None):
		self.filenames = Filenames(folder)

		self.label_order = label_order
		self.groups = groups

		self.anova_panel_plotter = AnovaPanelPlot(
			label_order_x = label_order,
			label_order_hue = groups
		)
		self.anova_plotter = AnovaPlotNested(
			label_order_x = label_order,
			label_order_hue = groups
		)

	def load(self):
		auc_statistics_table = pandas.read_csv(self.filenames.filename_table_auc_statistics, sep = "\t")

		return auc_statistics_table

	def run(self, ylimits:Tuple[int,int] = None):
		logger.debug(f"Saving the figures...")
		auc_statistics_table = self.load()
		is_nested = auc_statistics_table['condition'].nunique() != 1

		if not is_nested:
			self.anova_plotter.plot_single(auc_statistics_table, 'strain', 'auc_e', filename = self.filenames.filename_figure_anova_plot_main, ylims = ylimits)
		else:
			self.anova_plotter.plot(auc_statistics_table, 'condition', 'auc_e', 'strain', filename = self.filenames.filename_figure_anova_plot_main)
			self.anova_panel_plotter.anovaplotpanel(
				auc_statistics_table,
				x = 'condition', y = 'auc_e', hue = 'strain',
				filename = self.filenames.filename_figure_anova_plot_groups,
				control = 'RKS'
			)


def plot_qq(regression: linear_model.RegressionResultsWrapper, filename: Path):
	import matplotlib.pyplot as plt
	from statsmodels.graphics import gofplots
	gofplots.qqplot(regression.resid, fit = True, line = '45')
	plt.savefig(filename)


def plot_sigmas(sigmas: pandas.Series, filename: Path):
	other.plot_sigmas(sigmas, filename)


def plot_growthcurves(table: pandas.DataFrame, timeseries: pandas.DataFrame, folder_growthcurves: Path):
	growthcurve_plotter = PlotGrowthcurves(folder_growthcurves)
	growthcurve_plotter.plot_growthcurves(table, timeseries)


def save_auc_statistics_table(table: pandas.DataFrame, filename: Path):
	table.to_csv(filename, sep = "\t")


def save_table_tukey(tukey_results: Dict[str, Any], filename: Path, filename_json: Path):
	_temp = {key: utilities.tukey_to_json(value) for key, value in tukey_results.items()}

	filename_json.write_text(json.dumps(_temp, indent = 4, sort_keys = True))

	tables = list()
	for name, tukey_result in tukey_results.items():
		text = tukey_result.summary().as_csv()
		# Remove the header
		contents = StringIO("\n".join(text.split('\n')[1:]))
		df = pandas.read_csv(contents)
		df.columns = [i.strip() for i in df.columns]
		df['name'] = name
		try:
			df['pvalues'] = tukey_result.pvalues
		except AttributeError:
			message = "The version of TukeyHSD currently implemented in statsmodels does not have a `pvalues` attribute, so that will be missing from the result."
			logger.warning(message)
			df['pvalues'] = None
		tables.append(df)
	newtable = pandas.concat(tables)
	# Add duplicates for each series with reversed group1/group2 keys. This will make is easier to lookup a specific group or row.
	cleaner_tukey = CleanTukey()
	fulltable = cleaner_tukey.clean(newtable)
	fulltable.to_csv(filename, sep = '\t', index = False)
	return fulltable


def plot_tukey(tukey_results: Dict[str, Any], folder: Path, controls: Dict[str, str]):
	other.plot_tukey(tukey_results, folder, controls)


def save_tukey_matrix(table: pandas.DataFrame, folder_tukey, ext: str = '.tsv'):
	groups = table.groupby(by = "name")

	for name, group in groups:
		filename = folder_tukey / f"tukey.{name}{ext}"
		reverse_group = group.copy()
		reverse_group['group1'], reverse_group['group2'] = group['group2'].values, group['group1'].values
		df = group
		# df = pandas.concat([group, reverse_group])
		matrix = df.pivot(index = 'group1', columns = 'group2', values = 'meandiff').fillna(0)
		matrix.to_csv(filename, sep = "\t")


def save_anova(anova_table: pandas.DataFrame, filename: Path):
	anova_table.to_csv(filename, sep = '\t', index = True)


def save_maximum_growth(maximum_growth: pandas.Series, filename: Path):
	maximum_growth.sort_values().to_csv(filename, sep = '\t')


def save_table_growthcurve_models(table: pandas.DataFrame, filename: Path):
	table.to_csv(filename, sep = "\t")


def save_regression(regression: linear_model.RegressionResults, filename: Path):
	filename.write_text(str(regression.summary()))


def save_table_info(table_info: Dict[str, List[str]], filename: Path):
	with filename.open('w') as output:
		for key, strings in table_info.items():
			output.write(f"{key}\n")
			for s in strings:
				output.write(f"\t{s}\n")
