from pathlib import Path
import pytest
import pandas
from graphics import AnovaPlotNested
from loguru import logger
folder_data = Path(__file__).parent / "data"

@pytest.fixture
def auc_statistics_table()->pandas.DataFrame:
	filename = folder_data / "auc_statistics.tsv"
	return pandas.read_csv(filename, sep = "\t")

@pytest.fixture
def plotter()->AnovaPlotNested:
	obj = AnovaPlotNested(
		wildtype = 'WT',
		control = 'RKS',
		label_order_x = ["RKS", "Trp", "Fe3+", "Arg", "Asp", "Iso", "Lys", "Phe", "Met"],
		label_order_hue = ['WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA']
	)
	return obj
def get_mean(table:pandas.Series, variable:str, value:float)->float:
	t = table[table[variable]==value]
	return t['auc_e'].mean()
def test_calculate_mean_values_with_a_single_index(auc_statistics_table, plotter):

	conditions = auc_statistics_table['condition'].unique()
	logger.debug(conditions)
	groups = auc_statistics_table.groupby(by = ['condition', 'strain'])
	expected = groups.mean()['auc_e']

	for condition in conditions:
		table = auc_statistics_table[auc_statistics_table['condition'] == condition]
		expected_group = expected[condition]
		result = plotter._calculate_series_means(table, x = 'strain', y = 'auc_e', hue = 'strain')
		result.index = [i[0] for i in result.index]

		assert result.to_dict() == expected_group.to_dict()

def test_get_index_order_simple(auc_statistics_table, plotter):

	label_order_x = ['A', 'B', 'C']
	label_order_hue = ['1', '2', '3']

	expected = [
		('A', '1'), ('A', '2'), ('A', '3'),
		('B', '1'), ('B', '2'), ('B', '3'),
		('C', '1'), ('C', '2'), ('C', '3')
	]

	plotter.label_order_x = label_order_x
	plotter.label_order_hue = label_order_hue

	result = plotter._get_index_order()

	assert result == expected

	# Make sure the allowed_labels variables actually work.
	expected = [i for i in expected if 'C' not in i]
	result = plotter._get_index_order(allowed_x = ['A', 'B'])

	assert result == expected

def test_get_index_order_complex(plotter):

	allowed_x = ['Arg']
	allowed_hue = ['A244T', 'N274Y', 'N455K', 'P421L','WT', 'tRNA']

	expected = [('Arg', 'WT'), ('Arg', 'A244T'), ('Arg', 'N274Y'), ('Arg', 'N455K'), ('Arg', 'P421L'), ('Arg', 'tRNA')]
	result = plotter._get_index_order(allowed_x, allowed_hue)

	assert result == expected