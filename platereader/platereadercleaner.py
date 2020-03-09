import itertools
from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger
import re
FORMAT = "[strain].[condition].[plate].[replicate}]"
def is_strain(string: str) -> bool:
	""" Checks if the column has a valid name for a condition."""
	# Need to remove stuff like 'B2' or 'A13'
	pattern = "^[A-Z][0-9]{1,2}$"
	if re.search(pattern, string):
		return False
	elif string in {'Temp.', 'Cycle', 'Time'}:
		return False
	return True


class TableCleaner:
	def __init__(self):
		self.time_column_name_original = 'Time [s]'
		self.time_column_name = 'time'
		self.label_delimiter = '.'
		self.label_format = f"[strain]{self.label_delimiter}[condition]{self.label_delimiter}[plate]{self.label_delimiter}[replicate]"
		self.rename_columns_manually = False
		self.correct_rks = False

	@staticmethod
	def _group_columns_by_sample(columns: List[str]) -> Dict[str, List[str]]:
		groups = dict()
		for label in columns:
			parts = label.split(' ')

			if is_strain(parts[0]):
				strain = parts[0]
			else:
				strain = 'blank'
			if strain not in groups: groups[strain] = list()
			groups[strain].append(label)

		del groups['blank']  # Don't need this right now.

		return groups

	def _validate_column_label(self, label: str) -> str:
		""" Validates that the column follows the required format."""

		while True:
			label_parts = self._split_label(label)
			if label_parts is None and label != self.time_column_name:
				label = input(f"Rename '{label}' to follow the format '{self.label_format}':")

			else:
				break
		return label

	def _split_label(self, label: str) -> Optional[Tuple[str, str, str, str]]:
		try:
			strain, media, plate, replicate = label.strip().split(self.label_delimiter)
			return strain, media, plate, replicate
		except ValueError:
			return None

	@staticmethod
	def _get_sample_label_delimiter(label: str) -> str:
		if '.' in label:
			return '.'
		elif ' ' in label:
			return ' '
		else:
			message = f"Could not determine the delimiter for sample '{label}'"
			raise ValueError(message)

	def _rename_sample(self, label: str, plate: int) -> str:
		""" The original label may not conform perfectly to the expected specs."""
		try:
			delimiter = self._get_sample_label_delimiter(label)
			sample, *media, replicate = label.split(delimiter)
		except ValueError as exception:
			message = f"Could not extract the sample/replicate fields from '{label}'"
			logger.error(message)
			raise exception
		if self.correct_rks:
			# Ex. 'tRNA RKS Low PH 3'
			# logger.warning(f"Running patch to infer pH levels for '{label}'")
			lower = [i.lower() for i in media]
			if 'low' in lower:
				media = 'RKS-lowph'
			elif 'high' in lower:
				media = 'RKS-highph'
			else:
				media = 'RKS-neutral'
		else:
			media = delimiter.join(media)

		result = f"{sample}.{media.lower()}.{plate}.{replicate}"
		# logger.info(f"Renamed '{label}' to '{result}'")
		return result

	def _process_table_column_labels(self, table: pandas.DataFrame, plate: int) -> pandas.DataFrame:
		sample_groups = self._group_columns_by_sample(table.columns)
		#logger.info(f"Removing columns that do not correspond to a specific sample.")
		# logger.info(f"Will only keep {sorted(sample_groups.keys())}")

		# get rid of blank columns.
		sample_columns = sorted(itertools.chain.from_iterable(sample_groups.values()))
		reduced_table = table[[self.time_column_name_original] + sample_columns]

		# Rename the columns to remove whitespace and to add the replicate ID (plate vs. technical)
		# Patch to rename columns manually.
		reduced_table.columns = [self.time_column_name] + [self._rename_sample(i, plate) for i in sample_columns]

		if self.rename_columns_manually:
			logger.warning("Running patch to rename columns manually.")
			newcolumns = list()
			for column in reduced_table.columns:
				newlabel = input(f"Rename '{column}': ").strip()
				if newlabel != '':
					newcolumns.append(newlabel)
				else:
					newcolumns.append(column)

		import utilities
		utilities.validate_labels(reduced_table.columns)
		return reduced_table

	def _remove_empty_columns(self, table: pandas.DataFrame) -> pandas.DataFrame:
		# Also remove blank columns
		rejected_columns = [i for i in table.columns if 'blank' in i.lower()]
		for column in table.columns:
			if column == self.time_column_name: continue
			series = table[column]
			total = series.isna().sum()

			try:
				if total == len(series):
					rejected_columns.append(column)
			except ValueError as exception:
				logger.error(f"Found duplicate columns for '{column}'")
				logger.error(f"\tThe series type is {type(series)}")
				raise exception

		table = table[[i for i in table.columns if i not in rejected_columns]]
		return table

	def remove_redundant_time_columns(self, table: pandas.DataFrame) -> pandas.DataFrame:
		time_columns = table[[self.time_column_name]]
		table = table[[i for i in table.columns if i != self.time_column_name]]
		single_time_column = time_columns.iloc[:, 0]
		table[self.time_column_name] = single_time_column
		# Mave the time column so it is the first column.
		table = table[[self.time_column_name] + [i for i in table.columns if i != self.time_column_name]]
		return table

	def clean_table(self, table: pandas.DataFrame, plate: int) -> pandas.DataFrame:
		# Need to clean up the tables.
		# Remove the clomns that were empty and rename the columns to include replicate type.
		table = self._process_table_column_labels(table, plate)
		# table = self.rename_samples(table, plate)
		# Remove any extra stuff placed after the table
		# Use the time column to figure out when the table stops.
		# blank_rows = table[self.time_column_name].isna()
		# blank_rows = blank_rows[blank_rows].index
		# The first blank value represents the end of the table.
		# table = table.iloc[0:blank_rows[0]]

		# Convert the 'time' column to minutes and round.
		table.loc[:,self.time_column_name] = table[self.time_column_name].apply(lambda s: round(s / 60))
		table.loc[:,self.time_column_name] = table[self.time_column_name].astype(int)

		table = self._remove_empty_columns(table)
		# table = self.remove_redundant_time_columns(table)

		return table
