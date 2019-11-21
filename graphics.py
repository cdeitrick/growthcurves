from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import pandas
import seaborn
plt.style.use('ggplot')
import equations

from loguru import logger
class PlotGrowthcurves:
	def __init__(self, folder: Path, time_limit: Optional[int] = None):
		self.folder = folder
		self.time_limit = time_limit

		self.labelx = 'time (minutes)'
		self.labely = 'population'
		self.color_time_limit = 'tab:red'

		self.figure_format = '.svg'

	def plot_growthcurves(self, timeseries_table: pandas.DataFrame, coefficient_table: pandas.DataFrame):
		"""
		Parameters
		----------
		timeseries_table: The table from the plate reader. Column corresponds to a sample, each row corresponds to a timepoint.
		coefficient_table: A table (indexed by sample) with the fit parameters for each sample.
		"""

		# Group each technical replicate into a single graph.
		plt.style.use('fivethirtyeight')
		groups = coefficient_table.groupby(by = ['strain', 'condition', 'plate'])
		for label, coefficient_group in groups:
			label_str = '.'.join(label)
			filename = self.folder / f"growthcurve.{label_str}.{self.figure_format}"
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


class AnovaPlot:
	def __init__(self, wildtype:str, control:str):
		# Save the wild-type label so that it can be plotted as the first strain.
		self.control_strain = wildtype
		self.control_condition = control
		# The value to plot against the categorical variables.
		self.column_value = "auc_l"
		# The column to use as the x-label value
		self.column_condition = 'condition'
		# The column to use as the 'hue' variable.
		self.column_sample = 'strain'

		# Set the x,y,and hue variables to match the seaborn API.
		self.x = self.column_sample
		self.y = self.column_value
		self.hue = self.column_sample

		self.label_y = 'Fitness (AUC)'
		self.label_x = 'Strain'
		self.label_axis_fontsize = 36
		self.label_ticks_fontsize = 24
		self.number_of_columns = 3


	def anovaplot(self, growthcurves: pandas.DataFrame, x: str, y: str, hue: str, ax: Optional[plt.Axes] = None, filename: Optional[Path] = None) -> plt.Axes:
		""" Plots the fitness of strains (area under the curve) against the condition.
			Parameters
			----------
			growthcurves: pandas.DataFrame
				Should be indexed by sample id (ex. 'A244T.lle.1.1') and have an 'auc_l' column.
				The sample is is generated from the specific conditions that sample was grown in:
				ex. A244T.lle.1.1
					* strain A244T
					* condition lle
					* plate 1
					* technical replicate 1
			x: The category to plot on the x-axis
			y: The column to plot on the y-axis. Usually the area under the curve.
			hue: The column to determine series color from. May be the same as `x` for panels.
			ax: Optional[plt.Axes
				A preexisting plt.Axes object to add plots to.
			filename: Optional[Path]
		"""
		# plt.style.use('ggplot')

		if ax is None:
			plt.close()
			fig, ax = plt.subplots(figsize = (12, 10))
		# Make sure the wild-type variable is plotted first.
		if x == self.column_condition:
			xorder = [self.control_condition] + sorted(i for i in growthcurves[x].unique() if i != self.control_condition)
		elif x == self.column_sample:
			xorder = [self.control_strain] + sorted(i for i in growthcurves[x].unique() if i != self.control_strain)
		else:
			message = f"The selected x-axis category does not correspond to an existing column in the input table."
			logger.warning(message)
			xorder = None
		seaborn.stripplot(
			data = growthcurves,
			x = x, y = y, hue = hue,
			ax = ax,
			palette = 'Set1',
			dodge = True,
			order = xorder
		)
		ax.legend = False
		ax.set_ylim(growthcurves[y].min(), growthcurves[y].max())
		if filename:
			ymin = growthcurves[y].min() - 10
			ymax = growthcurves[y].max() + 10
			ax = self.format_plot_main(ax, (ymin,ymax))
			filename_svg = filename.with_suffix('.svg')
			filename_png = filename.with_suffix('.png')
			plt.savefig(filename_png, dpi = 500)
			plt.savefig(filename_svg)
		return ax
	def format_plot_main(self, ax:plt.Axes, ylims:Tuple[int,int] = None)-> plt.Axes:
		ax.set_xlabel("Condition", fontsize = self.label_axis_fontsize)
		ax.set_ylabel("Fitness (AUC)", fontsize = self.label_axis_fontsize)
		ax.tick_params(axis = 'both', labelsize = self.label_ticks_fontsize)

		if ylims:
			ax.set_ylim(ylims[0], ylims[1])
		return ax


	def _arrange_table(self, table: pandas.DataFrame, by:str, control:str) -> pandas.DataFrame:
		""" Basically sorts the data frame and put the control samples as the first sample in the table."""
		table_control = table[table[by] == control]
		table_notcontrol = table[table[by] != control]
		return pandas.concat([table_control, table_notcontrol])

	def anovaplotmultiple(self, growthcurves: pandas.DataFrame, groupby: str, filename: Optional[Path] = None):
		""" Plots multiple anova plots in the same figure. Each plot will corespond to a
			single value given by `groupby`.
			Parameters
			---------
			growthcurves: pandas.DataFrame
				A dataframe with the samples and corresponding auc measurements.
			groupby: str
				The variable to separate into multiple subplots within the panel.
			filename: Optional[Path]
				Where to save the figure.
		"""
		plt.style.use('fivethirtyeight')
		# Need to get the min/max auc values so that each subplot is plotterd useing an identical y-axis.
		ymin = growthcurves[self.column_value].min() - 10
		ymax = growthcurves[self.column_value].max() + 10  # Add 10 so that the largest value isn't partially off-screen

		unique_groups = growthcurves[groupby].unique()
		number_of_unique_groups = len(unique_groups)

		# Plot three plots per row.
		number_of_rows = int(number_of_unique_groups / self.number_of_columns) + 1

		figure: plt.Figure = plt.figure(figsize = (12, 9))

		grid = plt.GridSpec(number_of_rows, self.number_of_columns, hspace = 1)
		growthcurves = growthcurves.sort_values(by = [self.column_condition, self.column_sample])
		growthcurves = self._arrange_table(growthcurves, by = self.column_condition, control = self.control_condition)

		# move the control condition to the first position
		groups = growthcurves.groupby(by = groupby, sort = False)

		for index, (group_name, group) in enumerate(groups):
			# Determine the x-y index for gridspec from the `index` value.
			row, column = divmod(index, self.number_of_columns)
			current_ax: plt.Axes = figure.add_subplot(grid[row, column])
			ax = self.anovaplot(group, x = self.column_sample, y = self.column_value, hue = self.column_sample, ax = current_ax)

			# Format the current axis.
			self.format_plot(ax, group_name, (ymin, ymax), is_first = column != 0)
		if filename:
			filename_svg = filename.with_suffix('.svg')
			filename_png = filename.with_suffix('.png')
			plt.savefig(filename_png, dpi = 500)
			plt.savefig(filename_svg)

	def format_plot(self, ax: plt.Axes, name: str, ylims: Tuple[int, int], is_first: bool) -> plt.Axes:
		""" Applys formatting parameters to the plot.
			Parameters
			----------
			ax: plt.Axes
				The ax object with the plot.
			name: str
				The name/title of the ax.
			ylims: Optional[Tuple[int,int]]
				The minimum and maximum values observed over the dataset. Used to adjust the axis bounds so that the min and max points aren't
				cut off. Also used to make the y-axis consistent across all subplots in the panel.
			is_first:bool
				Whether this plot is the first one in the row. Used to determine whether to keep the y-axis labels/ticks
		"""

		if ylims:
			ax.set_ylim(ylims[0], ylims[1])
		# Remove the legend for individual plots. Save the first legend to add it later in an empty section of the panel.
		ax.get_legend().set_visible(False)
		ax.tick_params(axis = 'x', rotation = 45)
		ax.set_title(name, fontsize = 30)
		ax.set_xlabel(self.label_x, fontsize = 20)
		ax.set_ylabel(self.label_y, fontsize = 20)

		# Since each row shares the same y-axis, disable the y-axis in all but the first plot in each column.
		if is_first:
			ax.yaxis.set_ticklabels([])
			ax.set_ylabel("")

		return ax


def plot_sigmas(sigmas: pandas.Series, filename: Path):
	seaborn.distplot(sigmas)
	plt.savefig(filename)


def plot_tukey(tukey_results: Dict[str, Any], folder: Path):
	for name, tukey_result in tukey_results.items():
		filename = folder / f"tukey.{name}.png"
		ax = tukey_result.plot_simultaneous()
		plt.savefig(str(filename))


def plot_qq(residuals: pandas.Series):
	pass
