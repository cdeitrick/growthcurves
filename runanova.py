from io import StringIO
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import pandas
import statsmodels.api as sm
from loguru import logger
from statsmodels.regression import linear_model
from statsmodels.stats.multicomp import MultiComparison

import graphics
import growthcurver
import utilities


def tukeyhsd(statistics_table: pandas.DataFrame):
	subjects = ['plate', 'strain', 'condition']
	tukey_results = dict()
	for subject in subjects:
		tukey_result = MultiComparison(statistics_table['auc_l'], statistics_table[subject]).tukeyhsd()
		tukey_results[subject] = tukey_result

	statistics_table['condition:strain'] = statistics_table['condition'] + "-" + statistics_table['strain']
	mc = MultiComparison(statistics_table['auc_l'], statistics_table['condition:strain'])
	tukey_results['condition_strain'] = mc.tukeyhsd()

	return tukey_results


def anova(table: pandas.DataFrame) -> Tuple[linear_model.RegressionResults, pandas.DataFrame]:
	"""
		Calculates ANOVA
	Parameters
	----------
	table

	Returns
	-------
	model:
		*.params: A pandas.Series object with the calculated coefficients
	anova:

	"""
	# auc_aov <- aov(auc_l ~ condition*strain + plate, data=d_stat)
	regression = linear_model.OLS.from_formula('auc_l ~ condition*strain + plate', data = table).fit()

	anova_table = sm.stats.anova_lm(regression, typ = 1)
	return regression, anova_table


def convert_to_matrix(data: pandas.Series):
	""" Converts a 1D series (indexed by the key pair of labels) into a distance matrix"""
	pass


class Filenames:
	""" Holds the filenames for important files/figures/etc. Also implements the methods to save tables/figures."""

	def __init__(self, folder: Path, plot_growthcurves: bool = True, time_limit: Optional[int] = None):
		self.table_format = '.tsv'
		self.figure_format = '.svg'
		folder = utilities.checkdir(folder)
		self.folder_data = utilities.checkdir(folder / "data")
		self.folder_figure = utilities.checkdir(folder / "figures")

		# Tables
		self.table_statistics = self.folder_data / ("statistics" + self.table_format)
		self.table_growthcurve = self.folder_data / ("growthcurves" + self.table_format)
		self.table_regression = self.folder_data / "regression.txt"
		self.table_anova = self.folder_data / ("anova" + self.table_format)
		self.table_coefficients = self.folder_data / ("coefficients" + self.table_format)
		self.folder_tukey = utilities.checkdir(self.folder_data / "tukey")
		self.table_tukey = self.folder_tukey / ("tukey" + self.table_format)

		# Figures
		self.folder_figures_growthcurves = utilities.checkdir(self.folder_figure / "growthcurves")
		self.folder_figures_tukey = utilities.checkdir(self.folder_figure / "tukey")
		self.figure_sigmas = self.folder_figure / ("sigmas" + self.figure_format)
		self.figure_qq = self.folder_figure / ("qq" + self.figure_format)
		self.figure_anova_plot_main = self.folder_figure / ("anovaplot.main" + self.figure_format)
		self.figure_anova_plot_groups = self.folder_figure / ("anovaplot.panel" + self.figure_format)

		# Plotters
		self.plot_growthcurves = plot_growthcurves
		self.growthcurve_plotter = graphics.PlotGrowthcurves(self.folder_figures_growthcurves, time_limit = time_limit)

	def save_growthcurves(self, table: pandas.DataFrame, growthcurves: pandas.DataFrame):
		growthcurves.to_csv(self.table_growthcurve, index = True, sep = '\t')
		if self.plot_growthcurves:
			self.growthcurve_plotter.plot_growthcurves(table, growthcurves)

	def generate_tukey_matrix(self, table: pandas.DataFrame):
		groups = table.groupby(by = "name")

		for name, group in groups:
			filename = self.folder_tukey / f"tukey.{name}{self.table_format}"
			reverse_group = group.copy()
			reverse_group['group1'], reverse_group['group2'] = group['group2'].values, group['group1'].values

			df = pandas.concat([group, reverse_group])
			matrix = df.pivot(index = 'group1', columns = 'group2', values = 'meandiff').fillna(0)
			matrix.to_csv(filename, sep = "\t")

	def save_tukey_results(self, tukey_results: Dict[str, Any]):
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
		newtable.to_csv(self.table_tukey, sep = '\t', index = False)
		self.generate_tukey_matrix(newtable)
		graphics.plot_tukey(tukey_results, self.folder_figures_tukey)

	def save_sigmas(self, sigmas: pandas.Series):
		graphics.plot_sigmas(sigmas, self.figure_sigmas)

	def save_statistics_table(self, statistics_table: pandas.DataFrame):
		statistics_table.to_csv(self.table_statistics, sep = '\t', index = False)

	def save_regression(self, regression: linear_model.RegressionResults):
		self.table_regression.write_text(str(regression.summary()))

		data = {
			'coefficients': regression.params,
			'pvalues':      regression.pvalues,
		}
		df = pandas.DataFrame(data)

		df.to_csv(self.table_coefficients, sep = '\t', index = False)

	def save_qq(self, regression: linear_model.RegressionResults):
		import matplotlib.pyplot as plt
		from statsmodels.graphics import gofplots
		fig = gofplots.qqplot(regression.resid, fit = True, line = '45')
		plt.savefig(self.figure_qq)

	def save_anova(self, anova_table: pandas.DataFrame):
		anova_table.to_csv(self.table_anova, sep = '\t', index = False)


class GrowthCurveAnalysis:
	def __init__(self, project_folder: Path, plot_growthcurves: bool = True, wildtype:str = 'WT', control:str = 'RKS', time_limit: Optional[int] = None):
		self.project_folder = project_folder
		self.time_limit = time_limit
		self.time_column = 'Time'

		self.filenames = Filenames(self.project_folder, plot_growthcurves, self.time_limit)


		self.anova_plotter = graphics.AnovaPlot(wildtype, control)

	def load_data(self, filename: Path) -> pandas.DataFrame:
		"""
			Load the plate table with absorption values and clean up the data.
		Parameters
		----------
		filename
		"""
		table = pandas.read_excel(filename)

		# need to determine the name of the time column
		if 'Time' in table.columns:
			self.time_column = 'Time'
		elif 'time' in table.columns:
			self.time_column = 'time'
		else:
			message = f"Could not determine which column held the timepoints ('Time' or 'time'). Got {table.columns}"
			raise ValueError(message)
		utilities.validate_labels(table.columns)
		table = table.set_index(self.time_column)

		table = utilities.normalize_table(table)

		# Correct some labels
		table.columns = [i.replace('Iso', 'lle').replace("A224T", "A244T").replace('Rks', 'RKS') for i in table.columns]

		return table

	def run(self, filename: Path):
		# TODO: Renamve the variables more consistenly
		# ex. timeeries_table is saved as growthcurve table.
		# Also, the anova plotter requires the auc table.
		timeseries_table = self.load_data(filename)
		sample_metadata_table = utilities.extract_sample_metadata(timeseries_table)
		growthcurve_table = growthcurver.summarize_growth(timeseries_table.T, time_limit = self.time_limit)
		sample_table = sample_metadata_table.merge(growthcurve_table, left_index = True, right_index = True)

		auc_statistics_table = sample_metadata_table.merge(growthcurve_table, left_index = True, right_index = True)
		regression, anova_result = anova(auc_statistics_table)
		tukey_results = tukeyhsd(auc_statistics_table)
		# Generate the full anova plot
		self.anova_plotter.anovaplot(auc_statistics_table, x = 'condition', y = 'auc_l',hue =  'strain', filename = self.filenames.figure_anova_plot_main)
		# Generate the grouped anova plot.
		self.anova_plotter.anovaplotmultiple(auc_statistics_table, 'condition', self.filenames.figure_anova_plot_groups)

		self.filenames.save_sigmas(growthcurve_table['sigma'])
		self.filenames.save_anova(anova_result)
		self.filenames.save_regression(regression)
		self.filenames.save_statistics_table(auc_statistics_table)
		self.filenames.save_tukey_results(tukey_results)
		self.filenames.save_qq(regression)
		self.filenames.save_growthcurves(timeseries_table, sample_table)


def main():
	args = create_parser()
	analysis = GrowthCurveAnalysis(args.output, args.timelimit, wildtype = args.wildtype, control = args.control)
	analysis.run(args.filename)


def create_parser():
	import argparse
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"filename",
		help = "A table with a 'time' column and where columns correspond to samples and conditions" \
			   "The columns should be formatted as [strain].[condition].[plate].[replicate]",
		type = Path
	)
	parser.add_argument(
		"--output",
		help = "The folder to save all of the output files. If not given, an output folder will be generated based on the input filename.",
		type = Path,
		default = None
	)
	parser.add_argument(
		"--timelimit",
		help = "Sets the maximum time (in minutes) to include in the analysis. The platereader reports inconsistent measurements after ~41Hours, so adjust accordingly.",
		type = int
	)

	parser.add_argument(
		"--control",
		help = "The label applied to the control condition",
		type = str
	)

	parser.add_argument(
		"--wildtype",
		help = "The label of the wildtype strain.",
		type = str,
		default = 'WT'
	)


	args = parser.parse_args()
	return args


if __name__ == "__main__":
	main()
