from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
import pandas
from loguru import logger

try:
	from .anovaplot import AnovaPlotNested
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	from anovaplot import AnovaPlot


class AnovaPanelPlot:
	def __init__(self, label_order_x: Optional[List[str]] = None, label_order_hue: Optional[List[str]] = None):
		# Save the wild-type label so that it can be plotted as the first strain.
		# The value to plot against the categorical variables.
		# Create the plotter for the individual tables
		self.plotter = AnovaPlotNested(label_order_x = label_order_x, label_order_hue = label_order_hue)

		self.label_y = 'Fitness (AUC)'
		self.label_x = 'Strain'
		self.label_axis_fontsize = 36
		self.label_ticks_fontsize = 24
		self.number_of_columns = 3

	@staticmethod
	def _arrange_table(table: pandas.DataFrame, by: str, control: str) -> pandas.DataFrame:
		""" Basically sorts the data frame and puts the control samples as the first sample in the table."""
		table_control = table[table[by] == control]
		table_notcontrol = table[table[by] != control]
		return pandas.concat([table_control, table_notcontrol])

	def format_plot_main(self, ax: plt.Axes, ylims: Tuple[float, float] = None) -> plt.Axes:
		ax.set_xlabel("Condition", fontsize = self.label_axis_fontsize)
		ax.set_ylabel("Fitness (AUC)", fontsize = self.label_axis_fontsize)
		ax.tick_params(axis = 'both', labelsize = self.label_ticks_fontsize)

		if ylims:
			ax.set_ylim(ylims[0], ylims[1])
		return ax

	def format_plot(self, ax: plt.Axes, name: str, ylims: Tuple[float, float], is_first: bool) -> plt.Axes:
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

		ax.tick_params(axis = 'x', rotation = 45)
		ax.set_title(name, fontsize = 30)
		ax.set_xlabel(self.label_x, fontsize = 20)
		ax.set_ylabel(self.label_y, fontsize = 20)

		# Since each row shares the same y-axis, disable the y-axis in all but the first plot in each column.
		if is_first:
			ax.yaxis.set_ticklabels([])
			ax.set_ylabel("")

		return ax

	def _set_yaxis_limits(self, series: pandas.DataFrame):
		""" Sets self.ylimits using the min/max o fthe series values"""
		y_min = series.min()
		y_max = series.max()

		y_limits = (y_min, y_max)
		self.plotter.y_value_limits = y_limits

	def _calculate_number_of_rows(self, number_of_unique_groups: int) -> int:
		""" Calculates how many rows are needed to plot all of the unique plots. This assumes that three columns will be used."""
		number_of_rows = int(number_of_unique_groups / self.number_of_columns)  # + 1
		if number_of_rows == 0:
			number_of_rows = 1 # In case there are less than 4 plots
		return number_of_rows



	def anovaplotpanel(self, auc_statistics_table: pandas.DataFrame, x: str, y: str, hue: str, control: str, filename: Optional[Path] = None):
		""" Plots multiple anova plots in the same figure. Each plot will corespond to a
			single value given by `groupby`.
			Parameters
			---------
			auc_statistics_table: pandas.DataFrame
				A dataframe with the samples and corresponding auc measurements.
			x,y,hue
				The x, y and group columns
			filename: Optional[Path]
				Where to save the figure.
		"""
		plt.style.use('fivethirtyeight')
		# Need to get the min/max auc values so that each subplot is plotterd useing an identical y-axis.
		ymin = auc_statistics_table[y].min() - 10
		ymax = auc_statistics_table[y].max() + 10  # Add 10 so that the largest value isn't partially off-screen

		unique_groups = auc_statistics_table[x].unique()
		number_of_unique_groups = len(unique_groups)

		# Plot three graphs per row.
		number_of_rows = self._calculate_number_of_rows(number_of_unique_groups)+1

		if number_of_rows < 6:
			figsize = (12, 12)
		else:
			figsize = (12, 24)
		figure: plt.Figure = plt.figure(figsize = figsize)

		grid = plt.GridSpec(number_of_rows, self.number_of_columns, hspace = 1)

		auc_statistics_table = auc_statistics_table.sort_values(by = [x, hue])
		# Can probably remove this.
		auc_statistics_table = self._arrange_table(auc_statistics_table, by = x, control = control)

		# move the control condition to the first position
		groups = auc_statistics_table.groupby(by = x, sort = False)

		for index, (group_name, group) in enumerate(groups):
			# Determine the x-y index for gridspec from the `index` value.
			row, column = divmod(index, self.number_of_columns)
			logger.trace(f"row = {row}, column = {column}, number_of_columns = {self.number_of_columns}, number_of_rows = {number_of_rows}")
			current_ax: plt.Axes = figure.add_subplot(grid[row, column])

			# Plot the anovaplot on the given ax.
			self.plotter.label_order_x = [group_name]

			current_ax = self.plotter.plot(
				table = group,
				x = x, y = y, hue = hue,
				ax = current_ax,
				title = group_name
			)
			# Format the current axis.
			self.format_plot(current_ax, group_name, (ymin, ymax), is_first = column != 0)
		handles, labels = plt.gca().get_legend_handles_labels()
		by_label = dict(zip(labels, handles))
		plt.legend(
			by_label.values(), by_label.keys(), loc = 'lower right',
			#bbox_to_anchor = (.75, .25, .5, .5),
			#ncol = auc_statistics_table[hue].nunique(),
			framealpha = 0
		)
		if filename:
			filename_svg = filename.with_suffix('.svg')
			filename_png = filename.with_suffix('.png')
			plt.savefig(filename_png, dpi = 500)
			plt.savefig(filename_svg)

	def save_figure(self, ax: plt.Axes, filename: Path, ylims: Tuple[int, int]) -> plt.Axes:
		""" Saves the current figure as both a png and svg file."""
		ymin = ylims[0] * 0.95
		ymax = ylims[1] * 1.05

		ax = self.format_plot_main(ax, (ymin, ymax))

		filename_svg = filename.with_suffix('.svg')
		filename_png = filename.with_suffix('.png')
		plt.savefig(filename_png, dpi = 500)
		plt.savefig(filename_svg)

		return ax


