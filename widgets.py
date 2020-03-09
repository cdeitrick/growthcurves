from pathlib import Path
import pandas
from typing import *
import itertools
def squareform(self)->pandas.DataFrame:
	""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation.
	"""
	keys = sorted(set(itertools.chain.from_iterable(self.pairwise_values.keys())))
	_square_map = dict()
	for left in keys:
		series = dict()
		for right in keys:
			value = self.get(left, right, 0)
			series[right] = value
		_square_map[left] = series
	return pandas.DataFrame(_square_map)

if __name__ == "__main__":
	pass