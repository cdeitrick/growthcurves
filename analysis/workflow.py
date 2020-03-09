from pathlib import Path
from typing import *

import pandas
from loguru import logger
from statsmodels.sandbox.stats.multicomp import TukeyHSDResults  # Used to add a typing annotation to tukeyhsd()

import analysis
import projectoutput
import utilities
from analysis import growthcurver
from projectpaths import Filenames

TRACE = True
if TRACE:
	logger.remove()  # Need to remove the default sink so that the logger doesn't print messages twice.
	import sys

	logger.add(sys.stderr, level = "TRACE")

pandas.set_option('mode.chained_assignment', None)
EXPECTED_FORMAT = "[strain].[consition].[plate].[replicate]"


def extract_labels_from_metadata(metadata: pandas.DataFrame, column: str, allowed_strains: List[str] = None, allowed_conditions: List[str] = None) -> \
		List[str]:
	"""
		Extracts sample ids based on a specific column in the sample metadata table.
	Parameters
	----------
	metadata: pandas.DataFrame
	allowed_strains, allowed_conditions: List[str]
		The allowed strains and conditions.
	"""
	if allowed_strains is None:
		allowed_strains = metadata['strain'].unique()
	if allowed_conditions is None:
		allowed_conditions = metadata['condition'].unique()
	approved: List[str] = list()
	for index, row in metadata.iterrows():
		is_allowed_strain = row['strain'] in allowed_strains
		is_allowed_condition = row['condition'] in allowed_conditions
		if is_allowed_strain and is_allowed_condition:
			approved.append(index)
	return approved


def _extract_columns_by_media(columns: List[str], group: List[str]) -> List[str]:
	result = [i for i in columns if i.split('.')[1] in group]
	return result


class GrowthCurveAnalysis:
	def __init__(self, treatments: List[str] = None, strains: List[str] = None, time_limit: Optional[int] = None):
		self.time_limit = time_limit
		self.time_column = 'Time'

		self.treatments = treatments
		self.strains = strains

	@staticmethod
	def _remove_extra_labels(table: pandas.DataFrame, allowed_labels: List[str] = None) -> pandas.DataFrame:
		""" Removes columns from the table which are not to be included in the analysis."""
		if allowed_labels is None:
			return table
		passed_columns = list()

		for column in table.columns:
			for i in allowed_labels:
				if i.lower() in column.lower():
					passed_columns.append(column)

		return table[passed_columns]

	def set_paths(self, folder: Path):
		self.filenames = Filenames(folder)

	def prepare_table(self, table: pandas.DataFrame) -> pandas.DataFrame:
		"""
			Load the plate table with absorption values and clean up the data.
		Parameters
		----------
		table
		"""
		# need to determine the name of the time column
		if 'Time' in table.columns:
			self.time_column = 'Time'
		elif 'time' in table.columns:
			self.time_column = 'time'
		else:
			message = f"Could not determine which column held the timepoints ('Time' or 'time'). Got {list(table.columns)}"
			raise ValueError(message)
		table = table.set_index(self.time_column)

		table = utilities.normalize_table(table)

		return table

	@staticmethod
	def convert_letter_case(table: pandas.DataFrame) -> pandas.DataFrame:
		labelmap = {
			'rks': 'RKS',
			'arg': 'Arg',
			'asp': 'Asp',
			'fe3': 'Fe3+',
			'iso': 'Iso',
			'lys': 'Lys',
			'met': 'Met',
			'phe': 'Phe',
			'trp': 'Trp'
		}

		table.loc[:, 'condition'] = table['condition'].apply(lambda s: labelmap.get(s, s))
		return table

	def generate_growthcurve_table(self, table: pandas.DataFrame) -> pandas.DataFrame:
		growthcurve_timeseries_table_original = self.prepare_table(table)
		growthcurve_timeseries_table = self._remove_extra_labels(growthcurve_timeseries_table_original, self.treatments)

		return growthcurve_timeseries_table

	def save_results_tables(self, auc_statistics_table: pandas.DataFrame, anovaresults: pandas.DataFrame,
			regression: Any, tukey_results: Dict[str, TukeyHSDResults]):
		# projectoutput.save_table_info(table_info, self.filenames.filename_table_info)
		# projectoutput.save_maximum_growth(timeseries_table.max(), self.filenames.filename_table_maximum_growth)
		projectoutput.save_auc_statistics_table(auc_statistics_table, self.filenames.filename_table_auc_statistics)
		projectoutput.save_anova(anovaresults, self.filenames.filename_table_anova)
		projectoutput.save_regression(regression, self.filenames.filename_table_regression_model)

		tukey_table = projectoutput.save_table_tukey(tukey_results, self.filenames.filename_table_tukey, self.filenames.filename_data_tukey)
		projectoutput.save_tukey_matrix(tukey_table, self.filenames.folder_tables_tukey)
		projectoutput.plot_tukey(tukey_results, self.filenames.folder_figures_tukey, controls = {})

	def info(self, columns: List[str]) -> Dict[str, List[str]]:
		unique_strains = set(i.split('.')[0] for i in columns)
		unique_conditions = set(i.split('.')[1] for i in columns)
		unique_plates = set(i.split('.')[2] for i in columns)
		unique_replicates = set(i.split('.')[3] for i in columns)
		unique_combinations = set(".".join(i.split('.')[:2]) for i in columns)

		result = {
			'strains':      sorted(unique_strains),
			'conditions':   sorted(unique_conditions),
			'combinations': sorted(unique_combinations),
			'plates':       sorted(unique_plates),
			'replicates':   sorted(unique_replicates)
		}
		return result

	def summarize_growth(self, table: pandas.DataFrame) -> pandas.DataFrame:
		growthcurve_timeseries_table = self.generate_growthcurve_table(table)

		logger.info("Summarizing growth...")
		growthcurve_model_table = growthcurver.summarize_growth(growthcurve_timeseries_table.T, time_limit = self.time_limit)

		return growthcurve_model_table

	def run(self, table: pandas.DataFrame, auc_column: str, project_folder: Path = None):
		self.filenames = Filenames(project_folder)

		growthcurve_model_table = self.summarize_growth(table)
		sample_metadata_table = utilities.extract_sample_metadata(growthcurve_model_table.index)

		logger.info("Calculating auc statistics...")
		auc_statistics_table = sample_metadata_table.merge(growthcurve_model_table, left_index = True, right_index = True)

		regression, anova_result = analysis.anovanested(auc_statistics_table, auc_column)

		# Need to fix the labels in the AUC statistics table so they correctly formatted for the figures.
		auc_statistics_table = self.convert_letter_case(auc_statistics_table)

		logger.info("Running tukey...")
		tukey_results = analysis.tukeyhsd(auc_statistics_table, auc_column)

		logger.info("Saving tables...")

		self.save_results_tables(
			auc_statistics_table = auc_statistics_table,
			anovaresults = anova_result,
			regression = regression,
			tukey_results = tukey_results,
		)
		figure_workflow = projectoutput.FigureWorkflow(project_folder, self.treatments, self.strains)

		figure_workflow.run(ylimits = (0, auc_statistics_table['auc_e'].max()))
