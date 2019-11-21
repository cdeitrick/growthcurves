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
