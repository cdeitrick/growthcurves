from pathlib import Path
from typing import List, Optional

import pandas
from loguru import logger
from platereader.platereadercleaner import TableCleaner

FORMAT = "[strain].[condition].[plate].[replicate}]"


class PlateReaderParser:
	def __init__(self):
		self.plate: Optional[int] = None
		self.header_value = 'Cycle Nr.'  # Indicates that the line contains header values.
		self.time_column_name_original = 'Time [s]'
		self.time_column_name = 'time'
		# Assume that additional tables within the same file
		# correspond to the same samples, but was restarted 24 hours after the previous table.
		self.time_offset = 24 * 3600

		self.allowed_strains = {'WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA'}
		self.allowed_media = {'Arg', 'Asp', 'Fe3+', 'Phe', 'RKS'}

		self.cleaner = TableCleaner()

	def scan_for_typos(self, columns:List[str])->pandas.DataFrame:
		""" Searched for typos in the column names."""

		# Get the delimiter used for the columns
		delimiter = '.' if '.' in columns[0] else ' '
		for index, label in enumerate(columns):
			if delimiter not in label: continue
			try:
				strain, media, plate, replicate = label.split(delimiter)
			except ValueError as exception:
				message = f"Column '{label}' is not formatted as '{FORMAT}'"
				raise ValueError(message)
			istypo = (strain not in self.allowed_strains) or (media not in self.allowed_media)
			if istypo:
				logger.warning(f"Found a typo in column {index}: '{label}'")


	def _clean_subtable(self, subtable: pandas.DataFrame, columns: List[str], offset: int) -> pandas.DataFrame:
		""" Cleans up the subtables from a plate reader output file prior to combining them together."""

		subtable.columns = columns
		# Correct the timepoints. Each table represents a 24 hour chunk of time.
		subtable[self.time_column_name_original] = subtable[self.time_column_name_original] + offset
		return subtable

	def _combine_tables(self, table_list: List[pandas.DataFrame]) -> pandas.DataFrame:
		""" Combine all of the tables extracted from the plate reader.
		"""

		if len(table_list) == 1:
			return table_list[0]
		# Need to modify additional tables to correct the timepoints and column names.
		# Assume that only the first column has the proper column labels.
		first_table_columns = table_list[0].columns
		modified_tables = list()
		for index, table in enumerate(table_list):
			offset = self.time_offset * index
			subtable = self._clean_subtable(table, first_table_columns, offset)

			modified_tables.append(subtable)

		new_table = pandas.concat(modified_tables)

		return new_table

	def _get_table_indicies(self, table: pandas.DataFrame) -> List[int]:
		""" Extracts all tables from a plate reader output file, ignoring any metadata."""
		table_indicies_start = list()
		# Iterate over the dataframe explicitly
		# Quick to implement, and not slow enough to cause a problem.
		for index, row in table.iterrows():
			# Use the index value to skip the metadata rows.
			if row[0] == self.header_value:
				table_indicies_start.append(index)

		return table_indicies_start

	def _remove_blank_lines(self, minor_table: pandas.DataFrame) -> pandas.DataFrame:
		blank_lines = minor_table[self.time_column_name_original].isna()
		try:
			first_blank_line = blank_lines[blank_lines].index[0]
			minor_table = minor_table.iloc[0:first_blank_line]
		except IndexError:
			pass
		return minor_table

	def _extract_tables(self, filename: Path, indicies: List[int]) -> List[pandas.DataFrame]:
		""" Extracts all tables from a plate reader output file."""
		table_list = list()

		for start_index in indicies:
			minor_table = pandas.read_excel(filename, skiprows = start_index + 1)  # `skiprows` is 1-based rather than 0-based like the index.
			# The first NAN value in any of the columns indicates the end of the table.
			nonblank_table = self._remove_blank_lines(minor_table)
			table_list.append(nonblank_table)
		return table_list

	def read_table(self, filename: Path, plate: int):
		"""
			Reads in a tils growth curve table. The excel tables have a bunch of
			metadata before the table actually starts. Need to remove these lines.
			The plate reader can only loop through 24 hour intervals. Longer eperiments will result in multiple tables that
			need to be combined together.
		"""
		#logger.info(f"Reading {filename}")
		garbage_table = pandas.read_excel(filename)
		# Go through the first column and find the index of `Cycle Nr.`, which indicates the header row.
		# There may also be multiple tables.
		table_indicies_start = self._get_table_indicies(garbage_table)

		table_list = self._extract_tables(filename, table_indicies_start)
		combined_table = self._combine_tables(table_list)

		# Clean up the table
		# Need to change the column names

		clean_table = self.cleaner.clean_table(combined_table, plate)
		self.scan_for_typos(clean_table.columns)

		return clean_table
