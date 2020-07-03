"""
	Basically just lists the strains/conditions in a platereader table to help get a handle on which strains
	or conditions are available and whether there are some duplicates.

"""

from pathlib import Path
from typing import *
import pandas
from platereader.plateparser import PlateReaderParser

def get_columns(table:pandas.DataFrame)->List[str]:
	""" Looks for the `columns` row within a table and returns the strain/condition fields."""

	table = pandas.read_excel()

def get_header_row(table:pandas.DataFrame)->int:
	""" Returns the index of the row that contains the fieldnames."""
	# This is the very first field in the column names.
	key_value = "Cycle Nr."

	# So, locate the `key_value` in the first (leftmost) column
	# Since the column names are not the first row in the source tables we have to locate it using the column values.
	# Also, since the column names added by pandas are based on the first row they will basically just be
	# placeholders, so use the `columns` attribute to figure out where the first row is.

	first_column_label = table.columns[0]
	first_column = table[first_column_label]

	# It turns out that there's a *lot* of shitespace in the first column and it's not formatted well (.strip() won't work),
	# so instead just check if the `key_value` is contained in the cell.
	first_column = first_column.apply(lambda s: key_value in s if isinstance(s, str) else False)
	# There should only be one `True` value.
	index_of_header_row = first_column[first_column].index[0]
	return index_of_header_row


def main():
	filename = "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves/2020-04-06-growthcurves/source_tables/2.7.20.Asp.Phe.Rep1.xlsx"
	df = pandas.read_excel(filename)
	result = get_header_row(df)

	print(result)

	folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves/2020-04-06-growthcurves/source_tables/")

	parser = PlateReaderParser()
	for index, filename in enumerate(folder.iterdir()):
		if filename.suffix != ".xlsx": continue
		if 'arg' not in filename.name.lower(): continue
		#print(f"Reading {filename.name}")
		table = parser.read_table(filename, index)
		for column in table.columns:
			if 'arg' not in column.lower(): continue
			total = table[column].sum()
			print(filename.name)
			print(f"\t'{column}'\t{total}")

def main2():
	filename = "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves/2020-07-03-growthcurves/2020-07-03-growthcurves.tsv"
	table = pandas.read_csv(filename, sep = "\t")
	total = table.sum()
	print(total.to_string())
	import collections

	t = [i.split('.')[1] for i in table.columns if '.' in i]
	c = dict(collections.Counter(t))
	from pprint import pprint
	for k, v in sorted(c.items(), key =lambda s: s[1]):
		print(k,v)
	for column in table.columns:
		if 'arg' not in column.lower(): continue
		print(column)

if __name__ == "__main__":
	main2()