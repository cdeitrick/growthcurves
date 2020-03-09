from pathlib import Path
from graphics import AnovaPlotNested
import pandas
import matplotlib.pyplot as plt

folder = Path(__file__).parent.parent / "tests" / "data"
def get_table()->pandas.DataFrame:
	filename = folder / "auc_statistics.tsv"
	table = pandas.read_csv(filename, sep = "\t")
	return table

def get_plotter()->AnovaPlotNested:
	p = AnovaPlotNested(wildtype = 'WT', control = 'RKS',
		label_order_x = ["RKS", "Trp", "Fe3+", "Arg", "Asp", "Iso", "Lys", "Phe", "Met"],
		label_order_hue = ['WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA']
	)
	return p


def try_creating_a_grouped_anovaplot():

	table = get_table()
	plotter = AnovaPlotNested('WT', 'RKS')

	plotter.plot(table, x = 'condition', y = 'auc_e', hue = 'strain')
	plt.show()

def try_creating_an_anova_when_there_is_only_one_condition():
	table = get_table()
	plotter = get_plotter()

	table = table[table['condition'] == 'Fe3+']

	plotter.plot(table, x = 'strain', y = 'auc_e', hue = 'strain')
	plt.show()


if __name__ == "__main__":
	try_creating_an_anova_when_there_is_only_one_condition()