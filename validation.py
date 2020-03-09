from pathlib import Path
from pathlib import Path
from typing import Optional, Tuple, Union

import pandas
from loguru import logger

EXPECTED_FORMAT = "[strain].[consition].[plate].[replicate]"
AUC_COLUMN = 'auc_e'


class ValidateTable:
	# makes sure the table is formatted correctly.
	def __init__(self):
		self.time_column_label = "time"

		self.label_delimiter = '.'
		self.label_format = f"[strain]{self.label_delimiter}[consition]{self.label_delimiter}[plate]{self.label_delimiter}[replicate]"

		# The minimum growth that is expected. Used to check for no-growth
		self.minimum_growth = 0.1

	@staticmethod
	def read_table(filename: Union[str, Path]) -> pandas.DataFrame:
		if filename.suffix == '.csv':
			table = pandas.read_csv(filename)
		elif filename.suffix == '.tsv':
			table = pandas.read_csv(filename, sep = '\t')
		elif filename.suffix == '.xlsx' or filename.suffix == '.xls':
			table = pandas.read_excel(filename)
		else:
			message = f"Cannot determine the filetype of '{filename}'"
			raise ValueError(message)
		return table

	def check_table(self, table: Union[Path, pandas.DataFrame]) -> pandas.DataFrame:
		if not isinstance(table, pandas.DataFrame):
			# Assume it is a Pathlike object
			table = self.read_table(table)

		# Patch to remove extra columns
		logger.warning(f"Executing patch to remove extra columns.")
		rejected_columns = ['time.1', 'time.2']
		table = table[[i for i in table.columns if i not in rejected_columns]]

		# Check whether there is a time column
		table = self._check_for_time_column(table)

		# Make sure the time column is in minutes/seconds/whatever
		table = self._check_time_units(table)

		# Make sure the column labels are correctly formatted.
		table = self._check_column_labels(table)

		# Check for columns showing little to no growth.
		self._check_maximum_growth(table)

		return table

	def _check_maximum_growth(self, table: pandas.DataFrame) -> None:
		for column in table.columns:
			maximum_growth = table[column].max()
			if maximum_growth < self.minimum_growth:
				logger.warning(f"The sample '{column}' showed no growth ({maximum_growth} < {self.minimum_growth})")

	def _check_column_labels(self, table: pandas.DataFrame) -> pandas.DataFrame:
		new_columns = list()
		for column in table.columns:
			new_column = self._validate_column_label(column)
			new_columns.append(new_column)
		table.columns = new_columns

		return table

	def _validate_column_label(self, label: str) -> str:
		""" Validates that the column follows the required format."""
		while True:
			label_parts = self._split_label(label)
			if label_parts is None and label != self.time_column_label:
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
	def _check_time_units(table: pandas.DataFrame) -> pandas.DataFrame:
		logger.debug("ValidateTable._check_time_units() is currently disabled")
		return table

	def _check_for_time_column(self, table: pandas.DataFrame) -> pandas.DataFrame:
		columns = table.columns

		# Check for the expected time label
		if self.time_column_label in columns:
			new_columns = columns
		elif 'Time' in columns:
			new_columns = columns
			index = new_columns.index(self.time_column_label.capitalize())
			new_columns[index] = self.time_column_label
		else:
			message = f"Could not identify a time column from {list(columns)}"
			raise ValueError(message)
		table.columns = new_columns

		return table
