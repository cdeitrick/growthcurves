growthcurves

# Input Table
These scripts take a table where columns correspond to samples and rows correspond to specific timepoints
The timepoint column should be named either 'time' or 'Time'.
Each column label should be formatted as [strain].[condition].[plate].[replicate].

# Usage
```bash
python runanova.py
    --output /home/cld100/storage/projects/tils/growthcurves-2019-10-14/TilsGrowthcurves 
    --timelimit 2400 
    --contol RKS
    --wildtype WT
    /home/cld100/storage/projects/tils/growthcurves-2019-10-14/TilSGC.xlsx
```
