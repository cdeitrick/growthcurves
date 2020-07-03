from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import pandas
import seaborn

import utilities

plt.style.use('ggplot')
from analysis import equations

from tqdm import tqdm


class PlotGrowthcurves:
	def __init__(self, folder: Path, time_limit: Optional[int] = None):
		self.folder = folder
		self.time_limit = time_limit

		self.labelx = 'time (minutes)'
		self.labely = 'population'
		self.color_time_limit = 'tab:red'

		self.figure_format = 'png'  # TODO: make commandline option

	def plot_growthcurves(self,coefficient_table: pandas.DataFrame, timeseries_table: pandas.DataFrame):
		"""
		Parameters
		----------
		timeseries_table: The table from the plate reader. Column corresponds to a sample, each row corresponds to a timepoint.
		coefficient_table: A table (indexed by sample) with the fit parameters for each sample.
		"""

		# Group each technical replicate into a single graph.
		plt.style.use('fivethirtyeight')
		groups = coefficient_table.groupby(by = ['strain', 'condition', 'plate'])
		for label, coefficient_group in tqdm(groups, total = len(groups)):
			label_str = '.'.join(label)
			for filetype in ['png', 'svg']:
				f = utilities.checkdir(self.folder / filetype)
				filename = f / f"growthcurve.{label_str}.{self.figure_format}"
				group_timeseries = timeseries_table.T.loc[list(coefficient_group.index)]
				self.plot_group(group_timeseries, coefficient_group, filename)

	def plot_group(self, timeseries: pandas.DataFrame, fit_data: pandas.DataFrame, filename: Path):
		colormap = seaborn.color_palette('Paired', len(timeseries) * 2)

		fig, ax = plt.subplots(figsize = (12, 10))
		index = 0
		for sample_id, sample_timeseries in timeseries.iterrows():
			xdata = sample_timeseries.index
			ydata_empirical = sample_timeseries.values
			coefficients = fit_data.loc[sample_id]
			ydata_fit = [equations.logistic_equation(t, coefficients['k'], coefficients['N'], coefficients['r']) for t in sample_timeseries.index]

			ax.scatter(xdata, ydata_empirical, color = colormap[index], label = sample_id)
			ax.plot(xdata, ydata_fit, color = colormap[index + 1])
			index = index + 2
		if self.time_limit:
			ax.axvline(self.time_limit, color = self.color_time_limit, linestyle = ':', label = "time cutoff")
		ax.legend()
		ax.set_title(filename.stem)
		ax.set_xlabel(self.labelx)
		ax.set_ylabel(self.labely)
		plt.savefig(filename)
		plt.close(fig)
		return ax
