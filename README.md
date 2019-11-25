# growthcurves
These scripts are meant to perform an anova analysis based on strains of bacteria growth in various types of media while correcting "box effects" that arise from using multiple plates.
## Input Table
These scripts take a table where columns correspond to samples and rows correspond to specific timepoints
The timepoint column should be named either 'time' or 'Time'.
Each column label should be formatted as `[strain].[condition].[plate].[replicate]`.

## Usage
```bash
python runanova.py
    --output [output folder]
    --timelimit 2400 
    --contol RKS
    --wildtype WT
    [input table]
```

## Output

### Tables

#### `anova.tsv`
Presents a summary of the nested anova statistics from the input dataset.

#### `populationmodel.tsv`
Summarizes the statistics behind the anova analysis
Will have the following columns:

- `sample`: The sample id as it appears in the input table. It will be formatted as `[strain].[condition].[plate].[replicate]`.
- `strain`, `condition`, `plate`, `replicate`: Extracted from the `sample` column.
- `N`, `k`, `r`: The coefficients of a fitted logistic growth curve for the timeseries growth measurements behind each sample.
- `sigma`: The uncertainty behind the fitted growth curve.
- `auc_l`: The area under the curve from a logistic equation fitted to the growthcurve measurements for each sample.
-  `auc_e`: The empirical area under the curve calculated from the growth curve data directly. 

#### `regression.txt`
Presents the regression statistics behind the anova model. This is essentially a constructed linear model with the estimated coefficients for each variable in the anova model.

#### `maximumgrowth.txt`
Contains the maximum growth value observed for each sample, sorted from lowest to highest. Useful when checking for samples with little to no growth.

### Figures

#### `anovaplot.main`
   A figure presenting a general summary of the nested anova. It shows how the AUC values for each sample correlated with the sample, condition, and plate categorical variables.

### Output File Structure
```
.
|---- data
|----|---- tukey
|----|----|---- tukey.condition_strain.tsv
|----|----|---- tukey.condition.tsv
|----|----|---- tukey.plate.tsv
|----|----|---- tukey.strain.tsv
|----|----|---- tukey.tsv
|----|---- anova.tsv
|----|---- populationmodel.tsv
|----|---- regression.txt
|----|---- maximumgrowth.txt
|---- figures
|----|---- anovaplot.main.(png|svg)
|----|---- anovaplot.panel.(png|svg)
|----|---- qq.svg
|----|---- sigmas.svg
|----|---- tukey
|----|----|---- tukey.condition_strain.png
|----|----|---- tukey.condition.png
|----|----|---- tukey.plate.png
|----|----|---- tukey.strain.png
|----|---- growthcurves
|----|----|---- png
|----|----|---- svg
````
