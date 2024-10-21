

import pandas as pd


sample_table =  pd.read_table("SampleTable.txt")

sample_table['Pairs'] = sample_table['Fraction'].str.replace(' |fraction', '') + "_" + sample_table['Time'].astype(str).str.zfill(4) + "_" + sample_table['Dataset']
sample_table['TimePoints'] = sample_table['Fraction'].str.replace(' |fraction', '') + "_" + sample_table['Genotype'] + "_" + sample_table['Dataset']

print(sample_table['TimePoints'].value_counts())
print(sample_table['Pairs'].value_counts())

#print(list(sample_table['LibraryID'][sample_table['Pairs'] == 'pullout_00.0_dataset16']))
