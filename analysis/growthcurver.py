from pathlib import Path
from typing import Optional, Tuple

import numpy
import pandas
from loguru import logger
from scipy.integrate import trapz
from scipy.optimize import curve_fit

from analysis import equations


def summarize_growth(table: pandas.DataFrame, time_limit: Optional[int] = None) -> pandas.DataFrame:
	"""
		Fits the growth values to a logistic function.
		Assumes that `table` is formatted so that each row is indexed by sample.
	"""
	if time_limit:
		table = table[[i for i in table.columns if i <= time_limit]]

	results = list()

	for sample_name, sample_data in table.iterrows():
		normalized_data = sample_data - sample_data.min()
		xdata = normalized_data.index.values
		ydata = normalized_data.values
		try:
			# Include an initial guess for the parameters
			# This helps the curve fit to find the correct parameters without failing.
			(k, N, r), pcov = curve_fit(equations.logistic_equation, xdata, ydata, p0 = [1, .001, .004])
		except RuntimeError as exception:
			logger.warning(f"Could not process '{sample_name}'")
			raise exception

		result = {
			'sample': sample_name,
			'k':      k,
			'N':      N,
			'r':      r,
			'auc_l':  calculate_area_under_curve_ideal(max(normalized_data.index), k, N, r),
			'auc_e':  calculate_area_under_curve_empirical(normalized_data),
			'sigma':  calculate_goodness_of_fit(normalized_data, k, N, r)
		}
		results.append(result)
	df = pandas.DataFrame(results).set_index('sample')
	return df

def calculate_goodness_of_fit(empirical_data: pandas.Series, k: float, N: float, r: float) -> float:
	total = 0
	rdf = len(empirical_data) - 3
	for xpoint, ypoint in empirical_data.items():
		ipoint = equations.logistic_equation(xpoint, k, N, r)
		residual = (ypoint - ipoint) ** 2 / rdf
		if numpy.isnan(residual): continue
		total += residual
	return numpy.sqrt(total)


def calculate_area_under_curve_empirical(data: pandas.Series) -> float:
	""" Calculates the area under the curve of the descrete dataset."""
	return trapz(data.values, data.index)


def calculate_area_under_curve_ideal(t: int, k: float, N: float, r: float) -> float:
	area_under_curve = equations.logistic_equation_integral(t, k, N, r) - equations.logistic_equation_integral(0, k, N, r)
	return area_under_curve


def load_from_file(filename: Path) -> Tuple[float, float, float]:
	table = pandas.read_csv(filename, sep = '\t')
	# The values will be identical for all rows.
	row = table.iloc[0]
	return row['k'], row['n0'], row['r']


if __name__ == "__main__":
	pass
