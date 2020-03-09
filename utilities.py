from pathlib import Path
from typing import *
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


def extract_sample_metadata(sample_names:List[str]) -> pandas.DataFrame:
	# Attempts to extract strain, replicate, sample_type, etc from the input table.
	# Assume that column names are formatted as '{strain}.{condition}.{plate}.{replicate}'
	metadata_table = [get_sample_metadata(label) for label in sample_names]
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

def showcolumns(table:pandas.DataFrame):
	for column in sorted(table.columns):
		print(column)
# Convert this into a typeddict when available.
class TukeyFields:
	confint: List[float]
	data: Iterable[float]
	df_total: int
	groups: Iterable[str]
	groupsunique: Iterable[str]
	meandiffs: Iterable[float]
	pvalues: Iterable[float]
	q_crit: float
	reject: Iterable[bool]
	reject2: Iterable[bool]
	std_pairs: List[float]
	variance: float

def tukey_to_json(result)->Dict[str,Any]:
	""" Converts TukeyHSDResults to a dictionary."""

	data = {
		'confint': list(tuple(i) for i in result.confint),
		'data': list(float(i) for i in result.data),
		'df_total': int(result.df_total),
		'groups': list(result.groups),
		'groupsunique': list(result.groupsunique),
		'meandiffs': list(float(i) for i in result.meandiffs),
		'pvalues': list(float(i) for i in result.pvalues),
		'q_crit': float(result.q_crit),
		'reject': list(bool(i) for i in result.reject),
		'reject2': list(bool(i) for i in result.reject2),
		'std_pairs': list(float(i) for i in result.std_pairs),
		'variance': float(result.variance)
	}
	return data

if __name__ == "__main__":
	import datetime
	current_date = datetime.datetime.now()
	print(current_date.time())
	print(str(current_date.date()))