from pathlib import Path
from typing import *

import matplotlib.pyplot as plt
import pandas
import seaborn
from loguru import logger
import matplotlib

plt.style.use('fivethirtyeight')
plt.rcParams['svg.fonttype'] = 'none'
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
matplotlib.rcParams.update(new_rc_params)

class GroupPlot:
	def __init__(self):
		self.label_order_hue = "WT,A244T,N274Y,N455K,P421L,tRNA-Ile2".split(',')
		self.label_order_x = "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp".split(',')
		self.indexby = ['condition', 'strain']

		group_name_1 = 'Complements WT Defect'
		group_name_2 = 'Does not Complement WT Defect'

		self.treatment_groups = {
			group_name_2: ['Fe3+', 'Met', 'Trp'],
			group_name_1: ['Arg', 'Asp', 'Ile', 'Phe'],
			'Lys':        ['Lys'],
			'RKS':        ['RKS']
		}

		self.indicies = {
			'RKS':        0,
			"Lys":        1,
			group_name_1: slice(2, 5),
			group_name_2: slice(5, 9)
		}

		self.treatment_group_order = ['RKS', 'Lys', group_name_1, group_name_2]

		# Parameters controlling how the data is plotted.
		self.scale = 0.75
		self.dodge_value = 0.7

		self.marker_point = 'o'  # The marker to use when plotting the individual AUC values
		self.marker_mean = '_'  # The marker to use when plotting the mean of each group.
		self.marker_point_size = 5  # How large should the plotted points be?
		self.marker_point_alpha = 0.5  # Controls the transparency of the plotted markers.

		# Parameters to adjust the appearance of the plots
		# self.color_background = "#f0f0f0"  # Slightly grey
		self.color_background = '#FFFFFF'
		# The colors need to be assigned manually because while the `order` keyword correctly orders the x-axis labels, the color palette
		# Does not seem to adjust to the new order, so the `meanplot` in particular has different colors than the `dotplot` despite
		# presumably being plotted in an identical order.
		self.color_palette_key = 'tab10'
		self.color_axis = '#000000'  # The color of the x and y axes.
		self.legend_facecolor = "#FFFFFF"
		# Parameters that affect the whole figure
		self.figure_title = ""  # "Condition-Strain"
		#self.figure_xlabel = "Treatment"
		self.figure_xlabel = "Growth medium (RKS + supplement)"
		self.figure_ylabel = 'Fitness (AUC)'
		self.label_axis_size = 24
		self.title_size = 36
		self.label_legend_title = "Genotype"

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

		# ordered_index = self._get_index_order(group_means_x_values, group_means_hue_values)
		ordered_index = result.index
		try:
			result = result[ordered_index]
		except KeyError as exception:
			logger.error(f"There was an error when indexing the group means:")
			logger.error(f"The means were grouped as {list(result.index)}")
			logger.error(f"Attempted to reindex with {list(ordered_index)}")
			raise exception
		return result

	def add_figure_axis(self, subplots: Dict[str, plt.Axes]) -> Dict[str, plt.Axes]:
		""" Adds the x and y axes to the figure. These are currently based on the four subplots in the figure."""

		def set_spine_properties(spine):
			spine.set_visible(True)
			spine.set_edgecolor(self.color_axis)
			spine.set_linewidth(1)

		def disable_spine(spine):
			spine.set_linewidth(0)

		for label, axes in subplots.items():
			if label != 'RKS':
				axes.yaxis.set_ticklabels([])
				axes.set_ylabel("")
				disable_spine(axes.spines['left'])
				disable_spine(axes.spines['right'])
				disable_spine(axes.spines['top'])

			else:
				set_spine_properties(axes.spines['left'])
				disable_spine(axes.spines['right'])
				disable_spine(axes.spines['top'])
			set_spine_properties(axes.spines['bottom'])
			axes.tick_params(axis = 'y', which = 'minor', left = True)

		return subplots

	def add_annotations(self, annotation_ax: plt.Axes) -> plt.Axes:
		annotation_ax.yaxis.set_visible(False)
		annotation_ax.set_xticks([])
		annotation_ax.xaxis.set_tick_params(which = 'both', labelcolor = '#FFFFFF')
		annotation_ax.patch.set_alpha(0)
		annotation_ax.set_facecolor(None)
		annotation_ax.set_alpha(1)
		annotation_ax.set_xlabel(self.figure_xlabel, fontdict = {'fontSize': self.label_axis_size}, labelpad = -50)
		annotation_ax.set_frame_on(False)

		return annotation_ax

	def dotplot(self, table: pandas.DataFrame, x: str, y: str, hue: str, ax: plt.Axes) -> plt.Axes:
		"""
			Plots all of the `y` values across all categories in `x`.
		"""
		order = [i for i in self.label_order_x if i in table['condition'].unique()]
		# For some reason giving the table and specified columns via `x`, `y`, and `hue` makes the plot show and extra label for 'auc_e
		ax = seaborn.stripplot(
			#data = table,
			#x = x, y = y, hue = hue,
			x = table[x].tolist(),
			y = table[y].tolist(),
			hue = table[hue].tolist(),
			ax = ax,
			size = self.marker_point_size,
			palette = self.color_palette_key,
			alpha = self.marker_point_alpha,  # Make the points a little transparent so the mean value markers are easily visible.
			dodge = self.dodge_value,  # Enable for the main plot, disable for the panel plot.
			zorder = 1,  # Make sure these are plotted under the mean value markers.
			hue_order = self.label_order_hue,
			order = order
		)
		return ax

	def meanplot(self, table: pandas.DataFrame, x: str, y: str, hue: str, ax: plt.Axes) -> plt.Axes:
		""" Adds mean values for the categorical variables on the x-axis."""
		# group means is a pandas.Series object where the x-labels form the index and the values correspond to the means of the series.
		group_means: pandas.Series = self._calculate_series_means(table, y)
		# Make sure the index is in the right order.
		group_means = group_means.reset_index()
		order = [i for i in self.label_order_x if i in group_means['condition'].unique()]

		# order = self._get_index_order()
		# group_means = group_means.reindex(order).reset_index()
		# Since the group means are a series, it may be better to iterate over it and draw a bar at the corresponding locations.

		ax = seaborn.pointplot(
			data = group_means,
			x = x, y = y, hue = hue,
			ax = ax,
			palette = self.color_palette_key,
			dodge = self.dodge_value,
			order = order,
			hue_order = self.label_order_hue,
			join = False,  # Don't draw lines between the points.
			markers = self.marker_mean,
			scale = self.scale,
			ci = None,
			zorder = 2  # Make sure the means are plotted above the points so they aren't obscured.
		)
		return ax

	def format_subplot(self, ax: plt.Axes, x: str, ylimits: Optional[Tuple[float, float]] = None) -> plt.Axes:
		""" Adds labels to each axis and modifies to colorscheme a bit."""

		ax.set_facecolor(self.color_background)
		#ax.set_xlabel(x)
		ax.set_xlabel(" ")
		ax.set_ylabel(self.figure_ylabel, fontdict = {'fontSize': self.label_axis_size}, labelpad = -5)
		if ylimits:
			ax.set_ylim(*ylimits)
		ax.get_legend().remove()

		return ax

	def add_dotplots(self, table: pandas.DataFrame, figure: plt.Figure) -> Dict[str, plt.Axes]:
		""" Adds the individual scatterplots to the figure.
		"""

		grid = plt.GridSpec(8, 8)  # , hspace = 1)
		ylimits = (0, table['auc_e'].max() + 100)
		figure_axes = dict()
		for label in self.treatment_group_order:
			logger.info(f"Adding '{label}' to the plot.")
			categories = self.treatment_groups[label]

			_t = table[table['condition'].isin(categories)]

			current_ax = figure.add_subplot(grid[:-1, self.indicies[label]])
			current_ax = self.dotplot(table = _t, x = 'condition', y = 'auc_e', hue = 'strain', ax = current_ax)
			current_ax = self.meanplot(table = _t, x = 'condition', y = 'auc_e', hue = 'strain', ax = current_ax)
			self.format_subplot(current_ax, label, ylimits = ylimits)
			figure_axes[label] = current_ax


		return figure_axes


	def add_legend(self, ax: plt.Axes):
		""" ax: The current axis. Use `not-OAA-derived` for now so the legend is plotted in the bottom right of the whole figure."""
		handles, labels = ax.get_legend_handles_labels()
		logger.debug(labels)
		by_label = dict(zip(labels, handles))
		logger.debug(by_label)
		plt.legend(
			by_label.values(), by_label.keys(), loc = 'lower right',
			framealpha = 1,
			facecolor = self.legend_facecolor,
			title = self.label_legend_title
		)

	def plot(self, table: pandas.DataFrame, filename: Path = None):
		# Sort the table
		# table = table.sort_values(by = ["condition", "strain", "plate", "replicate"], ascending = False)

		# The last row in the grib will be used to add a third label.
		logger.debug(table.columns)
		table['strain'] = table['strain'].apply(lambda s: s.replace('A224T', 'A244T'))
		table['sample'] = table['sample'].apply(lambda s: s.replace('A224T', 'A244T'))
		table['strain'] = table['strain'].apply(lambda s: s.replace('tRNA', 'tRNA-Ile2'))
		table['sample'] = table['sample'].apply(lambda s: s.replace('tRNA', 'tRNA-Ile2'))
		grid = plt.GridSpec(8, 8)  # , hspace = 1)
		figure: plt.Figure = plt.figure(figsize = (12, 10))

		# figure.set_facecolor(self.color_background)
		# figure.set_edgecolor(self.color_background)

		figure.suptitle(self.figure_title, size = self.title_size)

		figure_axes = self.add_dotplots(table, figure)

		# Add the legend to the last plot.
		#logger.warning('patch')
		#self.treatment_group_order = [i for i in self.treatment_group_order if i in figure_axes]
		self.add_legend(figure_axes[self.treatment_group_order[-1]])

		self.add_figure_axis(figure_axes)

		# Add the annotation axes
		annotation_ax: plt.Axes = figure.add_subplot(grid[-1, :])
		self.add_annotations(annotation_ax)

		if filename:
			logger.info(f"Saving as {filename}.")
			plt.savefig(filename.with_suffix('.png'), facecolor = self.color_background)
			plt.savefig(filename.with_suffix('.svg'), facecolor = self.color_background)
		else:
			plt.show()


if __name__ == "__main__":
	tils_folder = Path.home() / "storage" / "projects" / "tils"
	growthcurve_folder = tils_folder / "growthcurves" / "debuggrowthcurves"
	# table_filename = growthcurve_folder / "2020-03-03-growthcurves.tsv"
	# /media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/growthcurves/2020-02-26-growthcurves/consolidated-notgrouped/data/auc_statistics.tsv

	table_filename = tils_folder / "growthcurves" / "2020-02-26-growthcurves/consolidated-notgrouped" / "data" / "auc_statistics.tsv"
	filename = growthcurve_folder / "condition_strain.png"
	t = pandas.read_csv(table_filename, sep = "\t")

	app = GroupPlot()
	app.plot(t)#, filename)
