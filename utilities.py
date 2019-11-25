from pathlib import Path
from typing import Dict, Union, List
from collections import Counter
import pandas
from loguru import logger
def checkdir(path:Union[str,Path])->Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path
def get_sample_metadata(sample_label: str) -> Dict[str, str]:
	""" Extracts metadata contained in the sample name.
		Assume the sample name is formatted as '{strain}.{condition}.{plate}.{replicate}'
	"""
	try:
		strain, condition, plate, replicate = sample_label.split('.')
	except ValueError:
		message = f"Could not parse a sample label ('{sample_label}'). Labels should be formatted as [strain].[condition].[plate].[replicate]"
		raise ValueError(message)
	data = {
		'sample':    sample_label,
		'strain':    strain,
		'condition': condition,
		'plate':     f"plate{plate}",
		'replicate': replicate
	}
	return data


def extract_sample_metadata(table: pandas.DataFrame) -> pandas.DataFrame:
	# Attempts to extract strain, replicate, sample_type, etc from the input table.
	# Assume that column names are formatted as '{strain}.{condition}.{plate}.{replicate}'
	metadata_table = [get_sample_metadata(label) for label in table.columns]
	return pandas.DataFrame(metadata_table).set_index('sample')


def normalize_series(data: pandas.Series) -> pandas.Series:
	""" Subtracts the minimum measured value from all other values in the series."""
	return data - data.min()


def normalize_table(table: pandas.DataFrame) -> pandas.DataFrame:
	""" Normalizes each column in the table by the minimum measured value for that column"""
	for column in table.columns:
		table[column] = normalize_series(table[column])
	return table

def _showcount(strings:List[str], index = None):
	""" Shows how often a string shows up."""

	if index is None:
		array = strings
	else:
		array = [i.split('.')[index] for i in strings]
	counter = Counter(array)

	for key, value in sorted(counter.items()):
		print(key,"\t", value)

def validate_labels(labels:List[str], verbose:bool = False):
	""" Checks all of the sample labels to make sure they conform to the format [strain].[condition].[plate].[replicate]"""
	# The 'time' column should be ignored for now.
	labels = [i for i in labels if i.lower() != 'time']
	for label in labels:
		parts = label.split('.')
		# Check for any additional parts in the sample label.
		if len(parts) != 4:
			message = f"The label '{label}' does not conform to the expected format: [strain].[condition].[plate].[replicate]"
			logger.warning(message)
	# Check for typos
	if verbose:
		print("Found the following occurances of strains: ")
		_showcount(labels, 0)

		print("Found the following occurances of conditions: ")
		_showcount(labels, 1)