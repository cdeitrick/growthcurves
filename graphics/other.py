from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
import pandas
import seaborn
from loguru import logger


def plot_sigmas(sigmas: pandas.Series, filename: Path):
	plt.close()
	plt.cla()
	plt.clf()
	seaborn.distplot(sigmas)
	plt.savefig(filename)


def plot_tukey(tukey_results: Dict[str, Any], folder: Path, controls:Dict[str,str]) -> Optional[plt.Axes]:
	"""
		Plots the various tukey plots.
	Parameters
	----------
	tukey_results
	folder
	controls: Dict[str,str]
		Maps the current subject to the control for that subject.

	Returns
	-------

	"""

	for name, tukey_result in tukey_results.items():
		filename = folder / f"tukey.{name}.png"
		try:
			ax = tukey_result.plot_simultaneous(comparison_name = controls.get(name))
			plt.savefig(str(filename))
		except Exception as exception:
			logger.error(exception)
			ax = None

	return ax
