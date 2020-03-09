from typing import *

import pandas
from pathlib import Path
import matplotlib.pyplot as plt

def info(table_columns: List[str]) -> Dict[str, List[str]]:
	columns = [i for i in table_columns if 'time' not in i]
	unique_strains = set(i.split('.')[0] for i in columns)
	unique_conditions = set(i.split('.')[1] for i in columns)

	unique_combinations = set(".".join(i.split('.')[:2]) for i in columns)

	result = {
		'strains':      sorted(unique_strains),
		'conditions':   sorted(unique_conditions),
		'combinations': sorted(unique_combinations)
	}
	return result

def convert_index(index:int, total_columns:int)->Tuple[int,int]:
	""" Converts an index to the row/column indicies."""
	index_row, index_column = divmod(index, total_columns)

	return index_row, index_column-1
def main():
	tils_folder = Path.home() / "storage" / "projects" / "tils"
	project_folder = tils_folder / "growthcurves" / "2020-02-26-consolidatedgrowthcurves"

	fname = "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves/2020-02-26-growthcurves/2020-02-26-growthcurves.tsv"
	table = pandas.read_csv(fname, sep = "\t")

	left = "A224T.Arg.1.1"
	right = "P421L.RKS.1.1"

	left = table[left]
	right = table[right]
	fig, ax = plt.subplots()
	ax.scatter(left.values, right.values)
	ax.plot([0,max(left.values)], [0,max(right.values)])
	plt.show()


if __name__ == "__main__":
	main()
