from pathlib import Path
from typing import *

import matplotlib.pyplot as plt

import pandas
import seaborn
from loguru import logger
import matplotlib

class AnovaPlotNested:
	""" Generates an anova plot with three variables: `x`, `y`, and `hue`"""

	def __init__(self, y_value_limits: Tuple[float, float] = None, label_order_x: Optional[List[str]] = None,
			label_order_hue: Optional[List[str]] = None):
		self.indexby = ['condition', 'strain']
		# Set the order in which the x-axis variables should be plotted
		self.label_order_x: Optional[List[str]] = label_order_x
		self.label_order_hue: Optional[List[str]] = label_order_hue

		# Parameters controlling how the data is plotted.
		self.scale = 0.75
		self.dodge_value = 0.7
		self.y_value_limits = y_value_limits

		self.marker_point = 'o'  # The marker to use when plotting the individual AUC values
		self.marker_mean = '_'  # The marker to use when plotting the mean of each group.
		self.marker_point_size = 5  # How large should the plotted points be?
		self.marker_point_alpha = 0.40  # Controls the transparency of the plotted markers.

		# Parameters to adjust the appearance of the plots
		self.color_background = "#f0f0f0"  # Slightly grey
		self.color_axis = '#000000'
		# The colors need to be assigned manually because while the `order` keyword correctly orders the x-axis labels, the color palette
		# Does not seem to adjust to the new order, so the `meanplot` in particular has different colors than the `dotplot` despite
		# presumably being plotted in an identical order.
		self.color_palette_key = 'tab10'
		if self.label_order_hue:
			self.color_palette = self._get_color_palette(self.label_order_hue)  # Will be generated when the plot is first created.
		else:
			self.color_palette = None

	def _get_color_palette(self, labels: List[str]) -> Dict[str, str]:
		""" Makes sure both plots are consistently colored. Has to be done since the x-axis labels are plotted in a different order."""
		color_palette = seaborn.color_palette(self.color_palette_key, len(labels))
		color_palette = {label: color for label, color in zip(labels, color_palette)}
		return color_palette

	def _get_index_order(self, allowed_x: Iterable[str] = None, allowed_hue: Iterable[str] = None) -> Union[List[str], List[Tuple[str, str]]]:
		if allowed_x is None:
			allowed_x = self.label_order_x
		else:
			allowed_x = list(allowed_x)
		if allowed_hue is None:
			allowed_hue = self.label_order_hue
		else:
			allowed_hue = list(allowed_hue)
		if self.label_order_hue is not None:
			result = list()
			for i in self.label_order_x:
				if i not in allowed_x: continue
				for j in self.label_order_hue:
					if j not in allowed_hue: continue
					result.append((i, j))
		else:
			result = [(i, i) for i in self.label_order_x if i in allowed_x]

		return result

	def get_colors(self) -> List[str]:
		""" Returns the color corresponding to the label."""
		labels = self.label_order_hue if self.label_order_hue is not None else self.label_order_x
		c = [self.color_palette.get(label, '#333333') for label in labels]

		return c

	def _calculate_series_means(self, auc_statistics_table: pandas.DataFrame, y: str) -> pandas.Series:
		""" This calculates the mean value of each x-hue pair. This results in a multi-indexed pandas.Series object,
			which is reduced to a single-dimension-ed index when `x` == `hue`
			Parameters
			----------
			auc_statistics_table: pandas.DataFrame
				The timeseries table.
		"""

		groups = auc_statistics_table.groupby(by = self.indexby)
		result = groups.mean()[y]
		group_means_x_values = set(i[0] for i in result.index)
		group_means_hue_values = set(i[1] for i in result.index)
		# if group_means_x_values == group_means_hue_values:
		#	result.index = [i[0] for i in result.index]

		ordered_index = self._get_index_order(group_means_x_values, group_means_hue_values)

		try:
			result = result[ordered_index]
		except KeyError as exception:
			logger.error(f"There was an error when indexing the group means:")
			logger.error(f"The means were grouped as {list(result.index)}")
			logger.error(f"Attempted to reindex with {list(ordered_index)}")
			raise exception
		return result

	def dotplot(self, table: pandas.DataFrame, x: str, y: str, hue: str, ax: plt.Axes) -> plt.Axes:
		"""
			Plots all of the `y` values across all categories in `x`.
		"""
		ax = seaborn.stripplot(
			data = table,
			x = x, y = y, hue = hue,
			ax = ax,
			size = self.marker_point_size,
			# palette = self.get_colors(),
			palette = self.color_palette_key,
			alpha = self.marker_point_alpha,  # Make the points a little transparent so the mean value markers are easily visible.
			dodge = self.dodge_value,  # Enable for the main plot, disable for the panel plot.
			order = self.label_order_x,  # The x-values are categorical, so they can be plotted in a specific order.
			zorder = 1,  # Make sure these are plotted under the mean value markers.

			# space = []
		)
		return ax

	def meanplot(self, table: pandas.DataFrame, x: str, y: str, hue: str, ax: plt.Axes) -> plt.Axes:
		""" Adds mean values for the categorical variables on the x-axis."""
		# group means is a pandas.Series object where the x-labels form the index and the values correspond to the means of the series.
		group_means: pandas.Series = self._calculate_series_means(table, y)
		# Make sure the index is in the right order.

		order = self._get_index_order()
		group_means = group_means.reindex(order).reset_index()
		# Since the group means are a series, it may be better to iterate over it and draw a bar at the corresponding locations.

		ax = seaborn.pointplot(
			# x = group_means.index,
			# y = group_means.values,
			data = group_means,
			x = x, y = y, hue = hue,
			ax = ax,
			# palette = self.color_palette_key,
			palette = self.get_colors(),
			dodge = self.dodge_value,
			# order = self.x_label_order,
			join = False,  # Don't draw lines between the points.
			markers = self.marker_mean,
			scale = self.scale,
			ci = None,
			zorder = 2  # Make sure the means are plotted above the points so they aren't obscured.
		)
		return ax

	def formatplot(self, ax: plt.Axes, x: str, y: str, title: Optional[str] = None, ylimits: Optional[Tuple[float, float]] = None) -> plt.Axes:
		""" Adds labels to each axis and modifies to colorscheme a bit."""
		if title:
			ax.set_title(title)

		ax.set_facecolor(self.color_background)
		if self.label_order_x and (self.label_order_x == 1):
			# This is probably one of the plots in the panalplot. So some of the labels can be omitted.
			ax.get_legend().remove()
			ax.xaxis.set_visible(False)
		else:
			# plt.legend(loc = 'lower right')
			handles, labels = plt.gca().get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			plt.legend(by_label.values(), by_label.keys(), loc = 'lower right', framealpha = 0)
			ax.set_xlabel(x)
			ax.set_ylabel(y)
		if ylimits:
			ax.set_ylim(*ylimits)

		return ax

	@staticmethod
	def save_figure(ax: plt.Axes, filename: Path) -> plt.Axes:
		""" Saves the current figure as both a png and svg file."""

		filename_svg = filename.with_suffix('.ps')
		filename_png = filename.with_suffix('.png')
		plt.savefig(filename_png, dpi = 500)
		plt.savefig(filename_svg)

		return ax

	def _reorder_table(self, table: pandas.DataFrame, x: str, hue: str) -> pandas.DataFrame:
		table_values_x = table[x].unique()
		table_values_hue = table[hue].unique()
		groups = table.groupby(by = self.indexby)
		dfs = list()
		label_order = self._get_index_order(allowed_x = table[x].unique(), allowed_hue = table[hue].unique())
		logger.trace(f"_reorder_table(, x = {x}, hue = {hue})")
		logger.trace(f"\tgroups: {groups.groups.keys()}")
		logger.trace(f"\tLabel_order_x: {self.label_order_x}")
		logger.trace(f"\tLabel_order_hue: {self.label_order_hue}")
		logger.trace(f"\ttable_values_x: {table_values_x}")
		logger.trace(f"\ttable_values_hue: {table_values_hue}")
		logger.trace(f"\tLabels: {label_order}")
		logger.trace(label_order)
		for element in label_order:
			group = groups.get_group(element)
			dfs.append(group)
		return pandas.concat(dfs)
	def add_figure_axis(self, subplots: Dict[str, plt.Axes]) -> Dict[str, plt.Axes]:
		""" Adds the x and y axes to the figure. These are currently based on the four subplots in the figure."""

		def set_spine_properties(spine):
			spine.set_visible(True)
			spine.set_edgecolor(self.color_axis)
			spine.set_linewidth(1)

		def disable_spine(spine):
			spine.set_linewidth(0)

		for label, axes in subplots.items():

			set_spine_properties(axes.spines['left'])
			disable_spine(axes.spines['right'])
			disable_spine(axes.spines['top'])
			set_spine_properties(axes.spines['bottom'])
			axes.tick_params(axis = 'y', which = 'minor', left = True)

		return subplots
	def plot_single(self, table: pandas.DataFrame, x: str, y: str, ax: Optional[plt.Axes] = None, filename: Path = None, ylims:Tuple[int,int] = None):
		plt.style.use('fivethirtyeight')
		plt.xticks(rotation = 70)
		if ax is None:
			fig, ax = plt.subplots(figsize = (10, 10))
		ax.set_facecolor(self.color_background)

		ax = seaborn.stripplot(
			data = table,
			x = x, y = y,
			ax = ax,
			size = self.marker_point_size,
			# palette = self.get_colors(),
			palette = self.color_palette_key,
			alpha = self.marker_point_alpha,  # Make the points a little transparent so the mean value markers are easily visible.
			dodge = False,  # Enable for the main plot, disable for the panel plot.
			order = self.label_order_x,  # The x-values are categorical, so they can be plotted in a specific order.
			zorder = 1,  # Make sure these are plotted under the mean value markers.
			# space = []
		)

		means = table.reset_index().groupby(by = 'strain').mean()['auc_e']

		ax = seaborn.pointplot(
			x = means.index,
			y = means.values,
			#data = means,
			#x = x, y = y,
			ax = ax,
			# palette = self.color_palette_key,
			palette = self.get_colors(),
			dodge = self.dodge_value,
			# order = self.x_label_order,
			join = False,  # Don't draw lines between the points.
			markers = self.marker_mean,
			scale = self.scale,
			ci = None,
			zorder = 2  # Make sure the means are plotted above the points so they aren't obscured.
		)
		ax = self.formatplot(ax, x = 'strain', y = 'auc_e')
		ax.set_title(table['condition'].unique()[0])
		ax = self.add_figure_axis({x:ax})[x]
		if ylims:
			ax.set_ylim(0, 1400)
		else:
			ax.set_ylim(0, table[y].max())
		ax.set_ylabel('Fitness (AUC)')
		ax.set_xlabel('WT Defect')
		ax.get_legend().remove()
		if filename:
			# Calculate the mean auc for each condition-strain value
			self.save_figure(ax, filename)
		return ax

	def plot(self, table: pandas.DataFrame, x: str, y: str, hue: str, ax: Optional[plt.Axes] = None, filename: Optional[Path] = None,
			title: Optional[str] = None, ylimits: Optional[Tuple[float, float]] = None):
		""" Plots the fitness of strains (area under the curve) against the condition.
			Parameters
			----------
			table: pandas.DataFrame
				Should be indexed by sample id (ex. 'A244T.lle.1.1') and have an 'auc_l' column.
				The sample is generated from the specific conditions that sample was grown in:
				ex. A244T.lle.1.1
					* strain A244T
					* condition lle
					* plate 1
					* technical replicate 1
			x: The category to plot on the x-axis
			y: The column to plot on the y-axis. Usually the area under the curve (AKA `auc_l`)
			hue: The column to determine series color from. May be the same as `x` for panels.
			ax: Optional[plt.Axes
				A preexisting plt.Axes object to add plots to.
			filename: Optional[Path]
			title: Optional[str]
				An optional title to add to the plot.
			ylimits: Optional[Tuple[float,float]]
				Used to make sure multiple plots can be made with the same y-scale.
		"""
		plt.style.use('fivethirtyeight')
		if x == hue:
			# "Since `x` and `hue` refer to the same variable, insert a copy of `x` into the table with the name `temp`.
			# This will prevent an error when the MultiIndex object created when calculating the series means is reset, since
			# that will throw an error when trying to reinsert the `x` column twice since they would have the same name.
			table['temp'] = table[x]
			hue = 'temp'
		self.indexby = [x, hue]

		table = table.sort_values(by = self.indexby)

		table = self._reorder_table(table, x, hue)
		plt.xticks(rotation = 70)
		if ax is None:
			fig, ax = plt.subplots(figsize = (10, 10))
		ax.set_facecolor(self.color_background)
		# Make sure the wild-type variable is plotted first.
		scatterplot = self.dotplot(table, x, y, hue, ax)

		# Show the conditional means
		meanplot = self.meanplot(table, x, y, hue, ax = scatterplot)

		meanplot = self.formatplot(meanplot, x, y, title, ylimits = ylimits)

		if self.y_value_limits:
			meanplot.set_ylim(self.y_value_limits[0], self.y_value_limits[1])

		if filename:
			# Calculate the mean auc for each condition-strain value

			self.save_figure(meanplot, filename)

		return meanplot
