from typing import *

import pandas
import statsmodels.api as sm
from loguru import logger
from statsmodels.regression import linear_model
from statsmodels.sandbox.stats.multicomp import TukeyHSDResults  # Used to add a typing annotation to tukeyhsd()
from statsmodels.stats.multicomp import MultiComparison


def tukeyhsd(statistics_table: pandas.DataFrame, column: str) -> Dict[str, TukeyHSDResults]:
	"""
		Perfors tukey multiple-comparison statistics.
	Parameters
	----------
	statistics_table: A table with each subject as a separate column
	column: The column with the relevant values. Should be identical to the `y` variable used when generating figures.
	"""
	is_nested = statistics_table['condition'].nunique() != 1
	if is_nested:
		subjects = ['plate', 'strain', 'condition']
	else:
		subjects = ['plate', 'strain']
	tukey_results = dict()
	for subject in subjects:
		logger.debug(f"tukey subject: '{subject}'")
		logger.debug(f"tukey subject values: {statistics_table[subject].unique()}")
		# MultiComparison doesn't work when there are only two possible groups, so disable this if we only have 2 categories.
		number_of_unique_categories = statistics_table[subject].nunique()
		if number_of_unique_categories > 2:
			tukey_result = MultiComparison(statistics_table[column], statistics_table[subject]).tukeyhsd()
			tukey_results[subject] = tukey_result

	statistics_table['condition:strain'] = statistics_table['condition'] + "-" + statistics_table['strain']
	mc = MultiComparison(statistics_table[column], statistics_table['condition:strain'])
	tukey_results['condition_strain'] = mc.tukeyhsd()

	return tukey_results


def anovanested(table: pandas.DataFrame, column: str) -> Tuple[linear_model.RegressionResults, pandas.DataFrame]:
	"""
		Calculates ANOVA
	Parameters
	----------
	table: The table containing the AUC values
	column: str; default 'auc_l'
		The column to get the AUC values from.

	Returns
	-------
	model:
		*.params: A pandas.Series object with the calculated coefficients
	anova:

	"""
	# auc_aov <- aov(auc_l ~ condition*strain + plate, data=d_stat)
	is_nested = table['condition'].nunique() != 1
	if is_nested:
		equation = f'{column} ~ condition + plate'
	else:
		equation = f'{column} ~ plate'
	logger.info(f"The equation used for ANOVA is {equation}")
	regression = linear_model.OLS.from_formula(equation, data = table).fit()

	anova_table = sm.stats.anova_lm(regression, typ = 1)
	return regression, anova_table



def main():
	from pathlib import Path
	project_folder = Path.home() / "storage" / "projects" / "tils" / "growthcurves" / "2020-02-26-growthcurves"

	filename = project_folder / "consolidated-notgrouped" / "data" / "auc_statistics.tsv"

	auc_statistics_table = pandas.read_csv(filename, sep = "\t")
	auc_column = "auc_e"
	auc_statistics_table = auc_statistics_table[auc_statistics_table['strain'] == 'WT']
	regression, anova_result = anovanested(auc_statistics_table, auc_column)
	"""
						df        sum_sq       mean_sq           F         PR(>F)
	condition           8.0  1.432147e+07  1.790184e+06  315.441447  3.539097e-206
	strain              5.0  2.463528e+06  4.927056e+05   86.817753   3.729951e-68
	plate              11.0  3.808582e+06  3.462347e+05   61.008684   1.596707e-89
	condition:strain   40.0  1.494660e+06  3.736650e+04    6.584207   1.830858e-27
	Residual          583.0  3.308625e+06  5.675171e+03         NaN            NaN
	
	"""

	print(anova_result)

if __name__ == "__main__":
	main()
