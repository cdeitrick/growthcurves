from pathlib import Path
import pytest
import pandas
from typing import *
from loguru import logger
from analysis import grouptools
folder_data = Path(__file__).parent / "data"
treatments = ['Arg', 'Fe3+']
samples = ['A244T', 'N274Y', 'N455K']
@pytest.fixture
def metadata()->pandas.DataFrame:
	# The auc table already includes the sample data so there's no point in having a redundant table in the data directory.
	filename = folder_data / "example_auc_statistics.tsv"

	table = pandas.read_csv(filename, sep = "\t").set_index('sample')

	table = table[['strain', 'condition', 'plate', 'replicate']]
	table = table[table['condition'].isin(treatments)]
	table = table[table['strain'].isin(samples)]

	return table

@pytest.fixture
def timeseries()->pandas.DataFrame:
	filename = folder_data / "example_timeseries.tsv"
	df = pandas.read_csv(filename, sep = "\t").set_index('time')


	return df

@pytest.fixture
def sample_columns()->List[str]:
	sample_columns = [
		'A224T.RKS.1.1',
		'A224T.Met.1.1',
		'A224T.Lys.1.1',
		'A224T.Asp.1.1',

		'N455K.RKS.1.1',
		'N455K.Met.1.1',
		'N455K.Lys.1.1',
		'N455K.Asp.1.1',

		'WT.RKS.1.1',
		'WT.Met.1.1',
		'WT.Lys.1.1',
		'WT.Asp.1.1',

		'P421L.RKS.1.1',
		'P421L.Met.1.1',
		'P421L.Lys.1.1',
		'P421L.Asp.1.1',
	]
	return sample_columns

def test_extract_labels_from_columns(metadata):
	result = workflow.extract_labels_from_metadata(metadata, column = 'condition', allowed_terms = [])

	expected = [
		'A224T.arg.1.1', 'A224T.arg.1.2', 'A224T.arg.1.3', 'A224T.fe3.1.1',	'A224T.fe3.1.2', 'A224T.fe3.1.3',
		'N274Y.arg.1.1', 'N274Y.arg.1.2', 'N274Y.arg.1.3',  'N274Y.fe3.1.1', 'N274Y.fe3.1.2','N274Y.fe3.1.3',
		'N455K.arg.1.1', 'N455K.arg.1.2', 'N455K.arg.1.3', 'N455K.fe3.1.1', 'N455K.fe3.1.2', 'N455K.fe3.1.3'
	]

	assert list(result) == list(expected)

@pytest.mark.parametrize(
	"strains, conditions, expected",
	[
		(["WT"], None, ['WT.RKS.1.1', 'WT.Lys.1.1', 'WT.Met.1.1', 'WT.Asp.1.1']),
		(["A224T", "N455K"], ["Met", "Lys"], ["A224T.Met.1.1", "A224T.Lys.1.1", "N455K.Met.1.1", "N455K.Lys.1.1"]),
		(['WT'], ['Met'], ['WT.Met.1.1'])
	]
)
def test_extract_sample_labels(sample_columns,strains, conditions, expected):
	# allowed_strains = ['A224T', 'WT', 'N455K']
	# allowed_conditions = ['RKS', 'Lys', 'Met']
	# Only compare against the samples representing plate 1 and replicate 1 to help simplify the test.


	result = grouptools.extract_sample_labels(sample_columns, strains, conditions)
	# Sort since the order doesn't matter, but the test checks for it.
	assert sorted(result) == sorted(expected)



def test_separate_columns_by_group(timeseries, sample_columns):
	groups = {
		'OTH': ['Lys', 'Met'],
		'RKS': ['RKS']
	}
	table = timeseries[[i for i in timeseries.columns if i.endswith('.1.1')]]

	result = grouptools.separate_columns_by_group(table.columns, groups)

	result = {key:sorted(value) for key, value in result.items()}


	expected_columns = {
		'OTH': [
			'A224T.Lys.1.1', 'N455K.Lys.1.1', 'WT.Lys.1.1',# 'P421L.Lys.1.1',
			'A224T.Met.1.1', 'N455K.Met.1.1', 'WT.Met.1.1',# 'P421L.Met.1.1',
		],
		'RKS': [
			'A224T.RKS.1.1', 'N455K.RKS.1.1', 'WT.RKS.1.1'#, 'P421L.Lys.1.1'
		]
	}
	expected_columns = {k:sorted(v) for k,v in expected_columns.items()}

	assert result == expected_columns

def test_group_by_plate_and_replicate(sample_columns):

	expected = {
		('A224T', '1', '1'): ['A224T.Lys.1.1', 'A224T.Met.1.1', 'A224T.RKS.1.1', 'A224T.Asp.1.1'],
		('N455K', '1', '1'): ['N455K.Lys.1.1', 'N455K.Met.1.1', 'N455K.RKS.1.1', 'N455K.Asp.1.1'],
		('WT', '1', '1'): ['WT.Lys.1.1', 'WT.Met.1.1', 'WT.RKS.1.1', 'WT.Asp.1.1'],
		('P421L', '1', '1'): ['P421L.Lys.1.1', 'P421L.Met.1.1', 'P421L.RKS.1.1', 'P421L.Asp.1.1']
	}

	result = grouptools.group_by_plate_and_treatment(sample_columns)
	expected = {key: sorted(value) for key, value in expected.items()}
	result = {key:sorted(value) for key, value in result.items()}
	assert result == expected

def test_calculate_mean_series(timeseries):
	series_1 = pandas.Series([1,2,3,4,5], name = 'series1')
	series_2 = pandas.Series([1,2,3,4,10], name = 'series2')
	series_3 = pandas.Series([1,1,1,1,1], name = 'series3')

	df = pandas.DataFrame([series_1, series_2, series_3]).transpose() # Ned to transpose so that each series forms a column rather than a row.

	expected12 = [1,2,3,4,7.5]
	expected13 = [1,1.5,2, 2.5, 3]

	result12 = grouptools.calculate_group_mean(df, ['series1', 'series2'])
	result13 = grouptools.calculate_group_mean(df, ['series1', 'series3'])
	result2 = grouptools.calculate_group_mean(df,['series2'])

	assert result12.tolist() == expected12
	assert result13.tolist() == expected13
	assert result2.tolist() == series_2.tolist()


if __name__ == "__main__":
	pass