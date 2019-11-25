import itertools
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger
FORMAT = "[strain].[consition].[plate].[replicate}]"

class PlateReaderParser:
	def __init__(self):
		self.plate: Optional[int] = None
		self.header_value = 'Cycle Nr.'  # Indicates that the line contains header values.
		self.time_column_name_original = 'Time [s]'
		self.time_column_name = 'time'
		# Assume that additional tables within the same file
		# correspond to the same samples, but was restarted 24 hours after the previous table.
		self.time_offset = 24 * 3600

		self.strains = {'WT', 'N455K', 'P421L', 'N274Y', 'tRNA', 'A224T'}

		self.cleaner = TableCleaner()

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
		logger.info(f"Reading {filename}")
		garbage_table = pandas.read_excel(filename)
		# Go through the first column and find the index of `Cycle Nr.`, which indicates the header row.
		# There may also be multiple tables.
		table_indicies_start = self._get_table_indicies(garbage_table)

		table_list = self._extract_tables(filename, table_indicies_start)

		combined_table = self._combine_tables(table_list)

		# Clean up the table
		# Need to change the column names
		#clean_table = self.clean_table(combined_table, plate)
		clean_table = self.cleaner.clean_table(combined_table, plate)
		return clean_table


class TableCleaner:
	def __init__(self):
		self.time_column_name_original = 'Time [s]'
		self.time_column_name = 'time'
		self.label_delimiter = '.'
		self.label_format = f"[strain]{self.label_delimiter}[consition]{self.label_delimiter}[plate]{self.label_delimiter}[replicate]"
		self.rename_columns_manually = False

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
	def _rename_sample(label: str, plate: int) -> str:
		""" The original label may not conform perfectly to the expected specs."""

		sample, *media, replicate = label.split(' ')

		logger.warning(f"Running patch to fix typo in '{label}'")
		lower = [i.lower() for i in media]
		if 'low' in lower:
			media = 'RKS-lowph'
		elif 'high' in lower:
			media = 'RKS-highph'
		else:
			media = 'RKS-neutral'

		result = f"{sample}.{media.lower()}.{plate}.{replicate}"
		logger.info(f"Renamed '{label}' to '{result}'")
		return result

	def _process_table_column_labels(self, table: pandas.DataFrame, plate: int) -> pandas.DataFrame:
		sample_groups = self._group_columns_by_sample(table.columns)
		logger.info(f"Removing columns that do not correspond to a specific sample.")
		logger.info(f"Will only keep {sorted(sample_groups.keys())}")

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
		table[self.time_column_name] = table[self.time_column_name].apply(lambda s: round(s / 60))
		table[self.time_column_name] = table[self.time_column_name].astype(int)

		return table


def is_strain(string: str) -> bool:
	""" Checks if the column has a valid name for a condition."""
	# Need to remove stuff like 'B2' or 'A13'
	pattern = "^[A-Z][0-9]{1,2}$"
	if re.search(pattern, string):
		return False
	elif string in {'Temp.', 'Cycle', 'Time'}:
		return False
	return True


def main_20191014():
	folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves-2019-10-14/original_tables")
	parser = PlateReaderParser()
	lys1 = folder / "TilSGC.Std.Lys.Iso.190815.xlsx"
	lys2 = folder / "TilSGC.Std.Lys.Iso.190822.xlsx"
	lys3 = folder / "TilsGC.Std.Lys.Iso.190830.xlsx"

	met1 = folder / "TilSGC.Std.Met.Trp.190906.xlsx"
	met2 = folder / "TilSGC.Std.Met.Trp.190929.xlsx"
	met3 = folder / "TilSGC.Std.Met.Trp.191008.xlsx"

	df_lys1 = parser.read_table(lys1, 1).reset_index(drop = True)
	df_lys2 = parser.read_table(lys2, 2).reset_index(drop = True)
	df_lys3 = parser.read_table(lys3, 3).reset_index(drop = True)
	df_met1 = parser.read_table(met1, 4).reset_index(drop = True)
	df_met2 = parser.read_table(met2, 5).reset_index(drop = True)
	df_met3 = parser.read_table(met3, 6).reset_index(drop = True)
	dfs = [df_lys1, df_lys2, df_lys3, df_met1, df_met2, df_met3]
	# Only need the 'Time' column ffor the first table.
	df_lys2 = df_lys2.drop(columns = ['time'])
	df_lys3 = df_lys3.drop(columns = ['time'])
	df_met1 = df_met1.drop(columns = ['time'])
	df_met2 = df_met2.drop(columns = ['time'])
	df_met3 = df_met3.drop(columns = ['time'])

	combined_table = pandas.concat([df_lys1, df_lys2, df_lys3, df_met1, df_met2, df_met3], axis = 1)
	combined_table = combined_table[['time'] + sorted(combined_table.columns[1:])]

	combined_table = combined_table.dropna()
	output_filename = folder.parent / "TilSGC.xlsx"
	logger.info(f"Saving as {output_filename}")
	combined_table.to_excel(output_filename, index = False)


def main():
	# TODO: remove extra time columns
	# TODO: Add a way to specify the expected strains and test for typos.
	project_folder = Path.home() / "storage" / "projects" / "tils" / "2019-11-22-growthcurves"
	table_folder = project_folder / "tables"
	output_filename = project_folder / "TilSGC.pH.tsv"
	parser = PlateReaderParser()

	tables = [
		table_folder / "TilSGC.Std.High.Low.191115.xlsx",
		table_folder / "TilSGC.Std.High.Low.191118.xlsx",
		table_folder / "TilSGC.Std.High.Low.191122.xlsx"
	]

	parsed_tables = [parser.read_table(filename, index + 1).reset_index(drop = True) for index, filename in enumerate(tables)]
	combined_table = pandas.concat(parsed_tables, axis = 1)
	combined_table = combined_table.dropna()
	combined_table.to_csv(output_filename, sep = '\t', index = False)


if __name__ == "__main__":
	main()
