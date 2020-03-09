import itertools
import math
from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
import pandas
import seaborn
from loguru import logger
class CleanTable:
	def __init__(self):
		self.strains = "WT,A244T,N274Y,N455K,P421L,tRNA".split(',')
		self.treatments = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')

	def extract_field(self, label: str, field: str, default: Any = None) -> Optional[str]:
		""" Tries to extract the condition label, if possible.

			Parameters
			----------
			label: str
				The label to classify
			field: {'condition', 'strain'}
				The field to get.
			default: Any
				What to return if the field is not found.
		"""
		if field == 'condition':
			index = 0
			allowed_fields = self.treatments
		else:
			index = 1
			allowed_fields = self.strains
		if '-' in label:
			# Corresponds to a condition_strain item. The treatment is the first argument, the strain is the second.
			result = label.split('-')[index]
		else:
			# Check if any of the treatments are in the label.
			try:
				result = [i for i in allowed_fields if i in label][0]
			except IndexError:
				result = default
		return result


	def run(self, table: pandas.DataFrame) -> pandas.DataFrame:
		table = self.add_inverse_keys(table)
		return table


class Application:
	def __init__(self, table:pandas.DataFrame):
		self.strains = "WT,A244T,N274Y,N455K,P421L,tRNA".split(',')
		self.treatments = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')
		self.table = table

	def _get_pairwise_treatment(self, table:pandas.DataFrame, field:str)->Dict[Tuple[str,str], float]:
		data = dict()
		for index, row in table.iterrows():
			key_left = row['group1']
			key_right = row['group2']
			value = row[field]
			if field == 'reject':
				value = int(value)
			data[key_left, key_right] = value
			data[key_right, key_left] = value
		return data

	def get_treatment_table(self, table:pandas.DataFrame, treatment:Tuple[str,str])->pandas.DataFrame:
		a = table['group1'].str.contains(treatment[0])
		b = table['group2'].str.contains(treatment[1])
		treatment_table = table[a & b]
		return treatment_table
	def plotsingle(self, table:pandas.DataFrame = None):
		if table is None:
			table = self.table

		print(table)

		pairwise_values = dict()
		table = table.set_index(['group1', 'group2'])

		series = table.loc[('Ile-WT', 'Asp-WT')]
		print(series)


	def plot(self):
		field = 'meandiff'
		table = self.table
		# Debug stuff
		size = len(table['group1'].unique())
		grid = plt.GridSpec(size, size)  # , hspace = 1)
		figure: plt.Figure = plt.figure(figsize = (12, 10))


		treatments = [(i,j) for i in table['condition1'].values for j in table['condition1'].values]
		for index, treatment in enumerate(treatments):

			treatment_table = self.get_treatment_table(self.table, treatment)
			pairwise_treatment = self._get_pairwise_treatment(treatment_table, field)
			square = squareform(pairwise_treatment)

			square = square[[i for i in square.columns if i.startswith(treatment[0])]]
			square = square.loc[[i for i in square.index if i.startswith(treatment[1])]]
			index_row, index_column = convert_index(index, len(self.treatments))
			ax:plt.Axes = figure.add_subplot(grid[index_row, index_column])

			palette = seaborn.color_palette("tab10")
			cmap = seaborn.diverging_palette(220, 10, as_cmap = True)
			seaborn.heatmap(square, ax = ax)
			ax.set_title("-".join(treatment))
		plt.show()



def convert_index(index: int, total_columns: int) -> Tuple[int, int]:
	""" Converts an index to the row/column indicies.
		from 2D to 1D:
			index = x + (y * width)
		from 1D to 2D:
			x = index % width
			y = index / width
	"""
	index_row = int(index / total_columns)
	index_column = index % total_columns

	return index_row, index_column


def tukey_to_squareform(table: pandas.DataFrame, field: str, keys: Tuple[str, str] = None) -> pandas.DataFrame:
	""" converts the tukey table to a pairwise table with values based on 'field'.
		The `keys` parameter selects which columns end up forming the dictionary key.
	"""
	if keys is None:
		keys = ['group1', 'group2']
	data = dict()
	seen = set()
	for index, row in table.iterrows():
		left = row[keys[0]]
		right = row[keys[1]]
		if (left, right) in seen: continue
		else:
			seen.add((left, right))
		value = row[field]
		data[left, right] = value

	result = squareform(data)
	return result


def squareform(pairwise_values: Dict[Tuple[str, str], float], default = math.nan) -> pandas.DataFrame:
	""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation.
	"""
	# Make sure the pairwise values are bidirectional
	parwise_values = pairwise_values.copy()
	#for key, value in pairwise_values.items():
	#	pairwise_values[key] = value
	#	pairwise_values[key[::-1]] = value
	labels = sorted(set(itertools.chain.from_iterable(pairwise_values.keys())))
	_square_map = dict()
	for left in labels:
		series = dict()
		for right in labels:
			value = pairwise_values.get((left, right), default)
			series[right] = value
		_square_map[left] = series
	return pandas.DataFrame(_square_map)


def compare(table: pandas.DataFrame, by: str) -> pandas.DataFrame:
	table = table.set_index(['group1', 'group2'])

	values = table[by]

	if by == 'reject':
		values = values.astype(int)

	square = squareform(values)
	return square


def main():
	project_folder = Path.home() / "storage" / "projects" / "tils" / "growthcurves" / "2020-02-26-growthcurves"
	filename = project_folder / "consolidated-notgrouped" / "data" / "tukey" / "tukey.tsv"
	other_filename = filename.with_name("tukey.condition_strain.tsv")
	auc_filename = filename.parent.parent / "auc_statistics.tsv"
	strains = "WT,A244T,N274Y,N455K,P421L,tRNA".split(',')
	treatments = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')
	treatment = "Lys"
	colors = {
		'A224T': "#e41a1c",
		'N274Y': "#377eb8",
		'N455K': "#4daf4a",
		'P421L': "#984ea3",
		'tRNA':  "ff7f00",
		'WT':    "#ffff33"
	}


	statistic_table = pandas.read_csv(filename, sep = "\t")
	statistic_table = statistic_table[statistic_table['name'] == 'condition_strain']
	statistic_table.columns = [i.strip() for i in statistic_table]

	wt = statistic_table[(statistic_table['group1'].str.contains('WT')) & (statistic_table['group2'].str.contains('WT'))]

	app = Application(wt)
	app.plotsingle(wt)



if __name__ == "__main__":
	main()
