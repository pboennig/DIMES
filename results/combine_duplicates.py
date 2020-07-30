import glob 
import pandas as pd
import re
import csv

duplicates = sorted(glob.glob("aucs/*_[12].csv"))

for i in range(0, len(duplicates), 2):
    cell_line = re.search("[^_]*", duplicates[i]).group(0)
    fst = pd.read_csv(duplicates[i], index_col=0)
    snd = pd.read_csv(duplicates[i+1], index_col=0)
    avg_auc = (fst['auc'] + snd['auc']) / 2
    fst['auc'] = avg_auc
    fst.to_csv(cell_line + ".csv", quoting=csv.QUOTE_NONNUMERIC)
