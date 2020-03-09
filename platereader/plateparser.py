from pathlib import Path
from typing import *

import pandas
from loguru import logger
import argparse
from platereader.platereaderparser import PlateReaderParser

FORMAT = "[strain].[condition].[plate].[replicate}]"


class Application:
	""" Basically acts as the entrypoint to the platereader parser."""

	def __init__(self):
		self.parser = PlateReaderParser()
		self.allowed_strains = []
		self.allowed_media = []

		self.labelmap = {
			'arg': 'Arg',
			'asp': 'Asp',
			'fe3': 'Fe3+',
			'iso': 'Ile',
			'lys': 'Lys',
			'met': 'Met',
			'phe': 'Phe',
			'rks': 'RKS',
			'trp': 'Trp'
		}

	def parse_tables(self, tables: List[Path], starting_plate: int = 0) -> List[pandas.DataFrame]:
		parsed_tables = list()
		for index, filename in enumerate(tables):
			plate_number = index + 1 + starting_plate
			logger.info(f"Converting {filename} representing plate {plate_number} to the required format...")
			result = self.parser.read_table(filename, plate_number).reset_index(drop = True)
			parsed_tables.append(result)
		return parsed_tables

	@staticmethod
	def summarize_table(table: pandas.DataFrame):
		""" Prints a list of all unique parameters in the table. This is usefull when looking for typos.
			Table
			-----
			columns: Formatted as '{strain}.{media}.{plate}.{replicate}'
		"""
		t = list()
		for column in table.columns:
			if column == 'time': continue
			try:
				strain, media, plate, replicate = column.split('.')
			except ValueError:
				message = f"Cannot extract metadata from label '{column}'"
				raise ValueError(message)

			r = {'strain': strain, 'media': media, 'plate': plate, 'replicate': replicate}
			t.append(r)
		t = pandas.DataFrame(t)
		for column in t.columns:
			counts = t[column].value_counts()
			print(f"Unique values found for '{column}'")
			for k, v in sorted(counts.items()):
				print(f"\t{v}\t{k}")

	def format_fields(self, labels:List[str])->List[str]:
		""" Corrects each column label using self.labelmap."""
		cols = list()
		for column in labels:
			for key, value in self.labelmap.items():
				if key in column:
					column = column.replace(key, value)
			cols.append(column)
		return cols
	def run(self, project_folder: Path, table_folder:Path, project_name: str = None):
		""" Provides the entrypoint for the platereader parser. This code here will likely be modified based on the current run."""
		# TODO: Add a way to specify the expected strains and test for typos.
		if project_name is None:
			project_name = project_folder.name
		output_filename = project_folder / f"{project_name}.tsv"

		tables_other = list(table_folder.iterdir())

		parsed_tables = self.parse_tables(tables_other)

		combined_table = pandas.concat(parsed_tables, axis = 1)
		combined_table = self.parser.cleaner.remove_redundant_time_columns(combined_table)


		# Convert to float?
		for column in combined_table.columns:
			if column == 'time': continue
			combined_table[column] = combined_table[column].astype(float)

		# Convert the condition fields from lowercase to tilte case.

		combined_table.columns = self.format_fields(combined_table.columns)
		self.summarize_table(combined_table)
		check_for_duplicate_columns(combined_table)
		combined_table.to_csv(output_filename, sep = '\t', index = False)

def lists_are_equal(left:List[float], right:List[float])->bool:
	import math
	results = list()
	for left_value, right_value in zip(left, right):
		result = math.isclose(left_value, right_value)
		results.append(result)

	return all(results)

def check_for_duplicate_columns(table:pandas.DataFrame):

	columns = table.columns
	seen = set()
	for left in columns:
		series_left = table[left].tolist()
		for right in columns:
			if left == right or (left, right) in seen: continue
			seen.add((left, right))
			seen.add((right, left))
			series_right = table[right].tolist()

			result = lists_are_equal(series_left, series_right)
			if result:
				print(f"Duplicate: '{left}' & '{right}'")


def main():
	# Correct some typos in the labels.

	project_name = "2020-02-26-growthcurves"
	project_folder = Path.home() / "storage" / "projects" / "tils" / "growthcurves" / project_name
	table_folder = project_folder /"tables" / "original_tables"

	Application().run(project_folder, table_folder)

def create_parser(args:Optional[List[str]] = None)->argparse.Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"folder",
		help = "The folder with the platereader tables.",
		type = Path
	)

	parser.add_argument(
		"--name",
		help = "The name of the project. Used to name the table.",
		type = str,
		default = None
	)

	args = parser.parse_args(args)

	return args

if __name__ == "__main__":
	main()
