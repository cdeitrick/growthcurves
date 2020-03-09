from pathlib import Path
from graphics import AnovaPanelPlot
import pandas
import matplotlib.pyplot as plt

folder = Path(__file__).parent.parent / "tests" / "data"
def get_table()->pandas.DataFrame:
	filename = folder / "auc_statistics.tsv"
	table = pandas.read_csv(filename, sep = "\t")
	return table

def get_plotter()->AnovaPanelPlot:
	p = AnovaPanelPlot(wildtype = 'WT', control = 'RKS',
		label_order_x = ["RKS", "Trp", "Fe3+", "Arg", "Asp", "Iso", "Lys", "Phe", "Met"],
		label_order_hue = ['WT', 'A244T', 'N274Y', 'N455K', 'P421L', 'tRNA']
	)
	return p

def test_that_it_works():
	table = get_table()
	plotter = get_plotter()

	result = plotter.anovaplotpanel(table, x = 'condition', y = 'auc_e', hue = 'strain')
	plt.show()