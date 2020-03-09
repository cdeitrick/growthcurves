from pathlib import Path
from typing import List

import pandas
import pytest

from platereader import plateparser


@pytest.fixture
def filename() -> Path:
	return DATA_FOLDER / "TilSGC.Std.Lys.Iso.190822.xlsx"


@pytest.fixture
def plate() -> plateparser.PlateReaderParser:
	return plateparser.PlateReaderParser()


@pytest.fixture
def subtables(filename, plate) -> List[pandas.DataFrame]:
	return plate._extract_tables(filename, [52, 219])


DATA_FOLDER = Path(__file__).parent / "data"


def test_get_table_indicies(filename, plate):
	expected = [52, 219]
	table = pandas.read_excel(filename)
	result = plate._get_table_indicies(table)
	assert result == expected


def test_remove_blank_lines(filename, plate):
	indicies = [52, 219]
	expected_index = list(range(1, 145))

	first_table = pandas.read_excel(filename, skiprows = indicies[0] + 1)
	# Make sure the first table still has blank lines
	assert first_table[plate.header_value].tolist() != expected_index

	nonblank_table = plate._remove_blank_lines(first_table)

	assert nonblank_table[plate.header_value].astype(int).tolist() == expected_index


def test_extract_tables(filename, plate):
	table = pandas.read_excel(filename)
	indicies = plate._get_table_indicies(table)

	tables = plate._extract_tables(filename, indicies)

	assert len(tables) == 2

	assert tables[0][plate.header_value].astype(int).tolist() == list(range(1, 145))
	assert tables[1][plate.header_value].astype(int).tolist() == list(range(1, 145))


def test_clean_subtable(subtables, plate):
	subtable_first, subtable_second = subtables
	columns = subtable_first.columns
	clean_subtable_first = plate._clean_subtable(subtable_first, columns, 0)

	assert list(clean_subtable_first.columns) == list(columns)
	assert clean_subtable_first[plate.time_column_name_original].values[0] == 0

	clean_subtable_second = plate._clean_subtable(subtable_first, columns, plate.time_offset*1)

	assert list(clean_subtable_second.columns) == list(columns)
	assert clean_subtable_second[plate.time_column_name_original].values[0] == plate.time_offset


def test_combine_tables(filename, plate):
	garbage_table = pandas.read_excel(filename)
	table_indicies_start = plate._get_table_indicies(garbage_table)
	table_list = plate._extract_tables(filename, table_indicies_start)

	combined_table = plate._combine_tables(table_list)

def test_scan_for_typos(filename, plate):
	pass