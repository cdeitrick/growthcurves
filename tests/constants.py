from pathlib import Path
from typing import *
from loguru import logger
import pandas
from toolz import itertoolz

folder_data = Path(__file__).parent / "data"
allowed_strains = ['A224T', 'WT', 'N455K']
allowed_conditions = ['RKS', 'Lys', 'Met']

def generate_example_timeseries_table(filename: Path, conditions: List[str], strains: List[str], timelimit = 100):
	logger.debug(f"generate_example_timeseries({filename.name}, {conditions}, {strains}, {timelimit})")
	output_filename = folder_data / "example_timeseries.tsv"

	table = pandas.read_csv(filename, sep = "\t").set_index('time')

	logger.debug(f"The input table has {len(table.columns)} columns.")
	cols = list()

	for column in table.columns:
		strain, condition, plate, replicate = column.split('.')
		#logger.debug(f"{strain, condition}")
		if strain in strains and condition in conditions:
			cols.append(column)
	logger.debug(f"Reduced the number of columns to {len(cols)}")
	table = table[cols]
	logger.debug(f"The second table has {len(table.columns)} columns.")
	table = table.loc[[i for i in table.index if i <= timelimit]]

	# Change the plate, replicate fields
	groups = itertoolz.groupby(lambda s: ".".join(s.split('.')[:2]), table.columns)

	wanted = ['.1.1', '.1.2', '.1.3', '.2.1', '.2.2', '.2.3', '.3.1', '.3.2', '.3.3']
	columnmap = dict()
	for name, group in groups.items():
		if name in ['A224T.RKS', 'N455K.RKS', 'WT.RKS']:
			result = {col:col for col in group[:9]}
		else:
			result = dict()
			for item_wanted, item_actual in zip(wanted, group):
				s,c,p,r = item_actual.split('.')
				item = f"{s}.{c}" + item_wanted
				result[item_actual] = item
		columnmap.update(result)

	# First remove any columns not in `columnmap`
	table = table[sorted(columnmap.keys())]

	# Rename columns
	table.columns = [columnmap[i] for i in table.columns]

	table.to_csv(output_filename, sep = "\t")


def generate_example_metadata_table(filename: Path, conditions: List[str], strains: List[str]):
	output_filename = folder_data / "example_auc_statistics.tsv"

	table = pandas.read_csv(filename, sep = "\t").set_index('sample')


	condition_series = table['condition'].isin(conditions)
	strain_series = table['strain'].isin(strains)

	table = table[condition_series & strain_series]


	table.to_csv(output_filename, sep = "\t")


def main():
	filename_timeseries = Path(
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/2020-02-21-growthcurves/2020-02-18-growthcurves.edited.tsv")
	filename_metadata = Path(
		"/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/2020-02-21-growthcurves/2020-02-21T16_48_14/data/auc_statistics.tsv")

	generate_example_timeseries_table(filename_timeseries, conditions = allowed_conditions, strains = allowed_strains)
	#generate_example_metadata_table(filename_metadata, conditions = allowed_conditions, strains = allowed_strains)


if __name__ == "__main__":
	main()
