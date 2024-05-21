"""
version: 1.21; 29-Mar-2023
this script need to run on server
"""
#%% 00 import statement
import pandas as pd
import cmapPy
from cmapPy import pandasGEXpress
from cmapPy.pandasGEXpress import parse
import pandas

#%% 01 parse LINCS gtcx
LINCS_chemicalPert_RNAseq = cmapPy.pandasGEXpress.parse.parse( 
    'LINCS_DCIC_2021_ChemicalPert_PredictedRNAseq_ChDir_Sigs.gctx') #size: 63.4G
print('parse LINS gtcx')

#%% 02 get raw dataframe
df_ChemicalPert = LINCS_chemicalPert_RNAseq.data_df
col_ChemicalPert = LINCS_chemicalPert_RNAseq.col_metadata_df
row_ChemicalPert = LINCS_chemicalPert_RNAseq.row_metadata_df
print('get raw dataframe')
    
#%% 03 set subset condition: cell. Here we use endometrosis cell lines as example
cell = pd.read_csv('uterus_cell_.csv').iloc[:, 0].tolist()

if len(cell) == 0:
    subset_col = col_ChemicalPert
else:
    subset_col = col_ChemicalPert.loc[col_ChemicalPert['cell'].isin(cell)]
    
#%% 04 get the sub-dataframe
subset_col = subset_col.transpose()
subset_df = pd.concat([subset_col, df_ChemicalPert], axis=0, join='inner', ignore_index=False)
print('get sub-dataframe')

#%% 05 export subset, if needed
pd.DataFrame.to_csv(subset_df, "sub_df_uterus.csv")
print('sub-dataframe exported')
