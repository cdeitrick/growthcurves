"""
	Generates a heatmap of heatmaps. The main heatmap is organized based on treatment, while each minor heatmap shows how the
	individual strains relate to each other within that treatment. It is easier to compare any combination of treatment/strain against all other
	combinations. The heatmap is based on the `tukey.tsv` file that was generated when each treatment/strain were treated independently.
"""

from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import pandas
import seaborn
from loguru import logger
import itertools
import matplotlib.patches as mpatches
from tqdm import tqdm
import matplotlib
plt.rcParams['svg.fonttype'] = 'none'
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
matplotlib.rcParams.update(new_rc_params)



def set_spine_properties(spine):
	spine.set_visible(True)
	spine.set_edgecolor('#333333')
	spine.set_linewidth(1)


def disable_spine(spine):
	spine.set_linewidth(0)


class Heatmap:
	def __init__(self):
		self.column_left = 'group1'
		self.column_right = 'group2'

		self.column_values = 'reject'  # The values that should be plotted in the individual heatmaps.
		self.color_axis = '#333333'

		self.color_positive = '#67a9cf'
		self.color_neutral = '#f7f7f7'
		self.color_negative = '#ef8a62'

		self.vmin = None
		self.vmax = None

		self.major = ['RKS', 'Trp', 'Fe3+', 'Met', 'Ile', 'Lys', 'Arg', 'Asp', 'Phe']
		self.minor = ['WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA']

	def generate_heatmap_minor(self, matrix: pandas.DataFrame, ax: plt.Axes, treatment_pair: Tuple[str, str]) -> plt.Axes:

		wt_label = [i for i in matrix.columns if 'WT' in i]
		other_labels = [i for i in matrix.columns if 'WT' not in i]
		matrix = matrix[wt_label + other_labels]

		wt_label = [i for i in matrix.index if 'WT' in i]
		other_labels = [i for i in matrix.index if 'WT' not in i]

		matrix = matrix.reindex(wt_label+other_labels)
		xticklabels = [i.split('-')[1] for i in matrix.columns]
		yticklabels = [i.split('-')[1] for i in matrix.index]

		ax = seaborn.heatmap(
			matrix,
			square = True,
			cmap = "RdBu_r",
			#cmap = [self.color_negative, self.color_neutral, self.color_positive][::-1],
			ax = ax, cbar = False,
			xticklabels = xticklabels, yticklabels = yticklabels,
			vmin = self.vmin, center = 0, vmax = self.vmax

		)
		ax.set_xlabel(treatment_pair[1], fontdict = {'size': 24})
		ax.set_ylabel(treatment_pair[0], fontdict = {'size': 24})

		# ax.set_title(treatment_pair)

		return ax

	@staticmethod
	def extract_treatment_table(table: pandas.DataFrame, treatments: Tuple[str, str], part:int = 0):
		rows = list()
		for index, row in table.iterrows():
			left = row['group1'].split('-')[part]
			right = row['group2'].split('-')[part]

			# if left in treatments and right in treatments:
			if left == treatments[0] and right == treatments[1]:
				rows.append(row)

		df = pandas.DataFrame(rows)

		return df

	def to_matrix(self, table: pandas.DataFrame, value_column: str, default: Any = 0) -> pandas.DataFrame:
		""" Converts a table with 'group1' and 'group2' columns into a matrix consisting of values from 'column'."""
		matrix = table.pivot(values = value_column, index = 'group1', columns = 'group2')
		matrix = matrix.fillna(default)
		return matrix

	def configure_axes(self, current_ax: plt.Axes, ax_position_index: Tuple[int, int], treatment_combinations_length: int) -> plt.Axes:
		# Configure the top axis
		if ax_position_index[0] == 0:
			set_spine_properties(current_ax.spines['top'])
			current_ax.xaxis.set_ticks_position('top')
			current_ax.xaxis.set_tick_params('both', rotation = 90)
			# current_ax.set_xlabel('X LABEL')
			current_ax.xaxis.set_label_position('top')
		elif ax_position_index[0] == treatment_combinations_length - 1:
			set_spine_properties(current_ax.spines['bottom'])
		else:
			current_ax.xaxis.set_visible(False)

		if ax_position_index[1] == 0:
			set_spine_properties(current_ax.spines['left'])
		elif ax_position_index[1] == treatment_combinations_length - 1:
			set_spine_properties(current_ax.spines['right'])
			current_ax.yaxis.set_ticks_position('right')
			current_ax.yaxis.set_tick_params('both', rotation = 0)
			current_ax.yaxis.set_label_position('right')
		else:
			current_ax.yaxis.set_visible(False)


		if ax_position_index[1] % 3 == 2:
			set_spine_properties(current_ax.spines['right'])
		elif ax_position_index[1] % 3 == 0:
			set_spine_properties(current_ax.spines['left'])

		if ax_position_index[0] % 3 == 2:
			set_spine_properties(current_ax.spines['bottom'])
		elif ax_position_index[0] % 3 == 0:
			set_spine_properties(current_ax.spines['top'])
		return current_ax

	def generate_masked_matrix(self, matrix_reject: pandas.DataFrame, matrix_values: pandas.DataFrame) -> pandas.DataFrame:
		matrix_values = matrix_values.copy(deep = True)
		masked_matrix = matrix_values.mask(~matrix_reject, 0)

		# Make sure half of the squareform matrix represents the opposite values.
		reduced_combinations = list()

		for index, row in masked_matrix.iterrows():
			values = [(index, col) for col in row[index:].index]
			reduced_combinations += values

		for i, c in reduced_combinations:
			value = masked_matrix.at[i,c]
			masked_matrix.at[i,c] = -value

		return masked_matrix

	def run(self, table: pandas.DataFrame, filename: Path):

		self.vmin = table['meandiff'].min()
		self.vmax = table['meandiff'].max()

		treatments = ['RKS', 'Trp', 'Fe3+', 'Met', 'Ile', 'Lys', 'Arg', 'Asp', 'Phe']
		treatments = ['WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA']
		treatment_combinations = [(i, j) for i in treatments for j in treatments]

		treatment_combinations_length = len(treatments)
		grid = plt.GridSpec(treatment_combinations_length, treatment_combinations_length)
		figure: plt.Figure = plt.figure(figsize = (20, 20))
		figure.suptitle("Comparison of fitness (AUC) by strain and treatment", size = 42)
		by_strain = int(True)
		for index, treatment_pair in enumerate(treatment_combinations):

			treatment_table = self.extract_treatment_table(table, treatment_pair, part = by_strain)
			print(treatment_pair)
			print(treatment_table)
			matrix_reject = self.to_matrix(treatment_table, 'reject', False)
			matrix_diff = self.to_matrix(treatment_table, 'meandiff', 0)

			masked_matrix = self.generate_masked_matrix(matrix_reject, matrix_diff)

			ax_position_index = get_position_index(index, treatment_combinations_length)
			current_ax = figure.add_subplot(grid[ax_position_index[0], ax_position_index[1]])
			current_ax = self.generate_heatmap_minor(masked_matrix, current_ax, treatment_pair)

			current_ax = self.configure_axes(current_ax, ax_position_index, treatment_combinations_length)
			if treatment_pair[0] == 'Trp' and treatment_pair[1] == 'Trp':
				print(masked_matrix.to_string())
		negative_patch = mpatches.Patch(facecolor = self.color_negative, label = 'X < Y', edgecolor = '#333333')
		neutral_patch = mpatches.Patch(facecolor = self.color_neutral, label = 'X = Y', edgecolor = '#333333')
		positive_patch = mpatches.Patch(facecolor = self.color_positive, label = 'X > Y', edgecolor = '#333333')
		plt.legend(
			handles = [negative_patch, neutral_patch, positive_patch],
			ncol = 3, loc = 'lower center',
			fontsize = 24,
			bbox_to_anchor = (-4.3, -1.5),
			frameon = False
		)


		if filename:
			plt.savefig(filename)
			plt.savefig(filename.with_suffix('.pdf'))
		plt.show()


def get_position_index(index: int, columns: int) -> Tuple[int, int]:
	index_column = index % columns
	index_row = int(index / columns)

	return index_row, index_column


def read_tukey_table(filename: Path, add_reverse:bool = False) -> pandas.DataFrame:
	t = pandas.read_csv(filename, sep = "\t")
	t = t[t['name'] == 'condition_strain']
	fulltable = t
	if add_reverse:
		fulltable['group1'], fulltable['group2'] = fulltable['group2'], fulltable['group1']
		fulltable['group1'] = fulltable['group1'].apply(lambda s: s.replace('A224T', 'A244T'))
	fulltable['group2'] = fulltable['group2'].apply(lambda s: s.replace('A224T', 'A244T'))
	fulltable['reject'] = fulltable['reject'].apply(lambda s: s.strip() == 'True')

	return fulltable


def main():
	tils_folder = Path.home() / "storage" / "projects" / "tils"
	project_folder = tils_folder / "heatmap"
	output_folder = project_folder / "figures"
	if not output_folder.exists():
		output_folder.mkdir()

	table_filename = project_folder / "tukey.tsv"
	t = read_tukey_table(table_filename, add_reverse = True)
	heatmap = Heatmap()
	heatmap.run(t, filename = None)#filename = output_folder / "heatmap.png")


if __name__ == "__main__":
	main()
