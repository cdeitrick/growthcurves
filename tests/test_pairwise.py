from pathlib import Path

import pandas
from analysis import pairwise
import pytest

@pytest.fixture
def tukey()->pandas.DataFrame:
	project_folder = Path.home() / "storage" / "projects" / "tils" / "growthcurves" / "2020-02-26-growthcurves"
	filename = project_folder / "consolidated-notgrouped" / "data" / "tukey" / "tukey.tsv"

	table = pandas.read_csv(filename, sep = "\t")

	return table

@pytest.fixture
def cleaner()->pairwise.CleanTable:
	return pairwise.CleanTable()

@pytest.mark.parametrize(
	"label, expected",
	[
		("Arg-A224T", "Arg"),
		("Ile-P421L", "Ile"),
		("Asp-WT", "Asp"),
		("Lys", "Lys"),
		("Fe3+", "Fe3+"),
		("A224T", None),
		("plate2", None)
	]
)
def test_extract_treatment_from_label(cleaner, label, expected):

	result = cleaner.extract_field(label, 'condition')
	assert result == expected
@pytest.mark.parametrize(
	"index, expected",
	[
		(0, (0,0)),
		(3, (1,0)),
		(10, (3,1)),
		(12, (4,0)),
		(7,  (2, 1))
	]
)
def test_convert_index(index, expected):
	"""
		(0,0)	(0,1)	(0,2)
		(1,0)	(1,1)	(1,2)
		(2,0)	(2,1)	(2,2)
		(3,0)	(3,1)	(3,2)
		(4,0)	(4,1)	(4,2)

	"""
	total_columns = 3


	result = pairwise.convert_index(index, total_columns)

	assert result == expected



if __name__ == "__main__":
	pass