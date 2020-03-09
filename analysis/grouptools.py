from typing import *

import pandas
from loguru import logger
from toolz import itertoolz

TRACE = False
if TRACE:
	logger.remove()  # Need to remove the default sink so that the logger doesn't print messages twice.
	import sys

	logger.add(sys.stderr, level = "TRACE")

pandas.set_option('mode.chained_assignment', None)
EXPECTED_FORMAT = "[strain].[consition].[plate].[replicate]"


def extract_sample_labels(labels: List[str], strains: List[str] = None, conditions: List[str] = None) -> List[str]:
	""" Extracts sample labels based on their strain or condition."""
	# logger.debug(f"allowed strains: {strains}")
	# logger.debug(f"allowed conditions: {conditions}")
	passed = list()
	for label in labels:
		strain, condition, plate, replicate = label.split('.')
		if strains:
			is_allowed_strain = strain in strains
		else:
			is_allowed_strain = True
		if conditions:
			is_allowed_condition = condition in conditions
		else:
			is_allowed_condition = True
		# logger.debug(f"{label}: {is_allowed_strain, is_allowed_condition}")
		if is_allowed_strain and is_allowed_condition:
			passed.append(label)
	return passed


def group_by_plate_and_treatment(columns: List[str]) -> Dict[Tuple[str, str, str], List[str]]:
	""" Combines all of the columns in `table` based on their `plate` and `replicate` values.
		Parameters
		----------
		columns: List[str]
			The columns to group
	"""

	def groupbykey(s) -> Tuple[str, str, str]:
		a, b, c, d = s.split('.')

		return a, c, d

	# Group the columns by plate, replicate, and strain while ignoring the condition.
	groups = itertoolz.groupby(groupbykey, columns)

	return groups


def separate_columns_by_group(columns: List[str], groups: Dict[str, List[str]]) -> Dict[str, List[str]]:
	""" Extracts the columns from `table` which are part of each group. Returns a dictionary mapping group names to the corresponding columns"""
	processed_group_tables = dict()
	for group_name, group_items in groups.items():
		logger.debug(f"Generating group {group_name}: {group_items}...")
		# Extract the sample ids based on the given criteria
		group_samples = extract_sample_labels(columns, conditions = group_items)
		logger.trace(f"Found {len(group_samples)} in group '{group_name}' that match this criteria.")

		processed_group_tables[group_name] = group_samples

	return processed_group_tables

def calculate_group_mean(table, columns:List[str])->pandas.Series:
	if len(columns) == 1:
		means = table[columns[0]]
	else:
		item = table[columns]
		means = item.mean(axis = 1)
	return means
def generate_grouped_series(table: pandas.DataFrame, groups: Dict[str, List[str]]) -> pandas.DataFrame:
	"""

	Parameters
	----------
	table: pandas.DataFrame
	groups: Dict[str,str]
		Maps groups of series under a common name.
	"""
	group_columns = separate_columns_by_group(table.columns, groups)

	# Each group is a dictionary mapping strain, plate, replicate to the corresponding columns.
	groups_by_plate_and_replicate:Dict[str,Dict[Tuple[str,str,str], List[str]]] = {key:group_by_plate_and_treatment(value) for key, value in group_columns.items()}

	# Go through each group and generate the mean table for each plate/replicate group.
	mean_table:Dict[str,pandas.Series] = dict()
	for group_name, group_data in groups_by_plate_and_replicate.items():
		for (strain, plate, replicate), columns in group_data.items():
			sample_name = f"{strain}.{group_name}.{plate}.{replicate}"
			mean_series = calculate_group_mean(table, columns)
			# `means` should be a pandas.Series object.
			mean_table[sample_name] = mean_series

	# Now combine the mean series for each group.
	df = pandas.DataFrame(mean_table)
	return df
