"""
	This file is only used as a reminder of how each table in the analysis is formatted
"""
from typing import Union

# Reminder of the format of each table.
class TableSchemaGrowthcurveModel:
	sample: str
	k:float
	N:float
	r:float
	auc_l:float
	auc_e:float
	sigma:float

class TableSchemaAucStatistics(TableSchemaGrowthcurveModel):
	# This table pairs the metadata for each sample with the fitted logistic curve for that sample.
	# This is identical to the Growthcurve model table, with extra columns with sample metadata
	strain: str
	condition:str
	plate:Union[str,int]
	replicate: Union[str,int]
	condition_strain:str # The actual field name is `condition:strain`

class TableSchemaAnova:
	# Contains the results of the ANOVA analysis
	df: int
	sum_sq:float
	mean_sq:float
	F:float
	PR:float # Actual name is PR(>F)

class TableSchemaTukey:
	group1: str
	group2: str
	meandiff: float
	p_adj:float # Actual name: p-adj
	lower: float
	upper: float
	reject:bool
	name: str # The category name. Ex. 'strain', 'plate'
	pvalues: float
