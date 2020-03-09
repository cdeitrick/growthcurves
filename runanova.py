from pathlib import Path
from typing import *

import pandas
from loguru import logger

import analysis
import utilities
from validation import ValidateTable

TRACE = False
if TRACE:
	logger.remove()  # Need to remove the default sink so that the logger doesn't print messages twice.
	import sys

	logger.add(sys.stderr, level = "TRACE")

pandas.set_option('mode.chained_assignment', None)
EXPECTED_FORMAT = "[strain].[consition].[plate].[replicate]"


def get_run_label() -> str:
	""" Generates the name of an output folder based on the current date and time."""
	import datetime
	current_date = datetime.datetime.now()
	date = str(current_date.date())
	time = current_date.time()
	# Rather than trying to delimit the time value, convert it to the total number of seconds after midnight.
	time_string = f"{time.hour}_{time.minute}_{time.second}"

	label = date + 'T' + time_string
	return label


def create_parser(args: List[str] = None):
	import argparse
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"filename",
		help = "A table with a 'time' column and where columns correspond to samples and conditions" \
			   "The columns should be formatted as [strain].[condition].[plate].[replicate]",
		type = Path
	)
	parser.add_argument(
		"--empirical",
		help = "Uses the empirical AUC values from the raw data rather than the integral of the fitted logistic curves.",
		action = "store_true"
	)
	parser.add_argument(
		"--output",
		help = "The folder to save all of the output files. If not given, an output folder will be generated based on the input filename.",
		type = Path,
		default = None
	)
	parser.add_argument(
		"--timelimit",
		help = "Sets the maximum time (in minutes) to include in the analysis. The platereader reports inconsistent measurements after ~41Hours, so adjust accordingly.",
		type = int
	)

	parser.add_argument(
		"--control",
		help = "The label applied to the control condition",
		type = str
	)

	parser.add_argument(
		"--wildtype",
		help = "The label of the wildtype strain.",
		type = str,
		default = 'WT'
	)
	parser.add_argument(
		"--treatments",
		help = "The order in which to plot the categorical variables on the x-axis. Sould be a comma-separated list.",
		type = str,
		default = None
	)
	parser.add_argument(
		"--strains",
		help = "The order in which to polt the subgroups within each `x` categorical variable. Should be a comma-separated list.",
		type = str,
		default = None
	)
	parser.add_argument(
		"--plot-growthcurves",
		help = "Whether to plot the measured values and fitted logistic equation for every sample. This may take a very long time.",
		action = "store_true",
		dest = "plotgrowthcurves"
	)
	if args:
		args = parser.parse_args(args)
	else:
		args = parser.parse_args()
	if args.treatments is not None:
		args.treatments = args.treatments.split(',')
	if args.strains is not None:
		args.strains = args.strains.split(',')
	return args


def main():
	project_folder = Path.home() / "storage" / "projects" / "tils" / "growthcurves" / "debuggrowthcurves"
	filename_table = project_folder / "2020-03-03-growthcurves.tsv"
	output_folder = utilities.checkdir(project_folder / "typetest")
	current_args = [
		'--output', str(output_folder),
		'--control', 'RKS',
		'--timelimit', '2400',
		# '--plot-growthcurves',
		'--empirical',
		'--treatments', "RKS,Lys,Arg,Asp,Fe3+,Ile,Met,Phe,Trp",
		'--strains', "WT,A244T,N274Y,N455K,P421L,tRNA",
		str(filename_table)
	]
	args = create_parser(current_args)

	validator = ValidateTable()
	table = validator.check_table(args.filename)
	analysis_workflow = analysis.GrowthCurveAnalysis(
		time_limit = args.timelimit,
		treatments = args.treatments,
		strains = args.strains

	)

	analysis_workflow.run(
		table,
		'auc_e' if args.empirical else 'auc_l',
		project_folder = project_folder / "debuggrowthcurves"
	)


if __name__ == "__main__":
	main()
