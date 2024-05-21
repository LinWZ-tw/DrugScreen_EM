"""
version: v2.2; 2023-10-13
author: lin, wei-zhi
"""
#%% 00
# 00 Import statement
from tqdm import tqdm
import pandas as pd
import numpy as np
import datetime
print ('00 modules imported')

#%% 01 identify soruce df # 01-1 set read csv func
def read_csv_with_progress(csv_file_path, chunk_size=100): #defaut chunk size is 100
    """
    Read a CSV file into a DataFrame with a progress bar.
    Parameters:
        csv_file_path (str): Path to the CSV file to be read.
        chunk_size (int, optional): The number of rows to read at a time. Defaults to 1000.
    Returns:
        pd.DataFrame: The DataFrame containing the data from the CSV file.
    """
    # Read the CSV file using tqdm.pandas() to get the progress bar
    tqdm.pandas()
    
    # Initialize the progress bar with an estimated total based on the chunk size
    with tqdm(total=None, unit=' lines', dynamic_ncols=True) as pbar:
        df_list = []
    
        # Use chunksize to read the CSV file in chunks
        for chunk in pd.read_csv(csv_file_path, chunksize=chunk_size):
            df_list.append(chunk)
            pbar.update(len(chunk))
            
    # Concatenate all chunks into a single DataFrame
    df = pd.concat(df_list, ignore_index=True)
    print ()
    return df

print ('01-1 func prepared')
#%% 01-2 Get LINCS df
""" source of df is generated from 1_parsing_LINCS """
#file_path_df = 'sub_df_uterus_test.csv' # this is a test file
file_path_df = 'sub_df_uterus.csv' # this is for real
print ('Set dataframe file path, loading...')
df = read_csv_with_progress(file_path_df).set_index('Unnamed: 0')
df_header = pd.DataFrame([df.columns], columns=df.columns)

# Insert header as a new row at the beginning
df = pd.concat([df_header, df]).rename(columns=df.iloc[1])

print ('01-2 Get LINCS dataframe')

#%% 01-3 Get DEG
""" difference gene expression from GEPIA2 in which the source data were from TCGA """
file_path_DEG = 'GEPIA2_DGA_UCEC_q=0.01_20230719.csv'
DEG = read_csv_with_progress(file_path_DEG)

print ('01-3 Get DEG')

#%% 02 filter gene by difference gene expression between disease tissue and normal tissue
""" set a threshold for q-value"""
q_threshold = 1e-4
# filter row that with adjp lower than q_threshold
DEG = DEG[DEG.adjp < q_threshold]
# get a list of gene in DEG
gene = DEG.iloc[:, 0].tolist()
gene_len = len(gene)
print ('02 Complete gene filtering')

#%% 03 Establish a dict
# filter row in df that with gene in DEG
df_col = df.iloc[:6]
df_filter = df.loc[df.index.isin(gene)]
df_filter = pd.concat([df_col,df_filter]).rename(columns=df.iloc[1]) # define cell line as header
cell_values = df.iloc[1].values # identify  cell values

# create a dict
dict_groups ={}

# group based on cell_values
for cell in set(cell_values):
    cell_mask = cell_values == cell
    cell_df = df_filter.loc[:, cell_mask]
    cell_df.columns = cell_df.iloc[0]

    dict_groups[cell] = cell_df

print ('03 Establish dict group by cell')

df = pd.DataFrame() # empty the unused df to release RAM (when nessesary)

#%% 04 calculate DEG[log2(fold change)] + df in dict_groups
# 04-1 Modify DEG to make index as the same as dict_groups
gene = df_filter.index[7:].tolist() # make sure gene start from 7
DEG_sub = DEG.loc[DEG['Gene Symbol'].isin(gene)].set_index('Gene Symbol')
DEG_value = pd.to_numeric(DEG_sub['Log2(Fold Change)'])

print ('04-1 Complete DEG modifying')


# 04-2 concate dict_groups and DEG_sub
dict_DEG = {}

for cell, cell_df in dict_groups.items(): 
    print ('processing', cell)
    dict_DEG[cell] = {}  # Create an empty dictionary for this cell
    temp_df = pd.DataFrame() # Create an empty df to store calculate result
  
    for column in cell_df.columns:
        #print (column)
        temp_DEG = cell_df[column].iloc[7:].astype(float)
        """ 
        Calculate the sum of 'gene expression difference by treatment' and 'DEG', using abs() 
        """
        temp_result = pd.DataFrame({f'{column}_result': (temp_DEG + DEG_value).abs()})
        temp_DEG = pd.DataFrame() # empty the unused df to release RAM
        
        temp_df = pd.concat([temp_df, temp_result], axis=1)

    # Store the result DataFrame in the dictionary
    dict_DEG[cell] = temp_df
    temp_df = pd.DataFrame() # empty the unused df to release RAM

print ('04-2 Calculate the sum of gene expression difference by treatment and DEG')


# 04-3 renamed the header of cell_df
for cell, cell_df in dict_DEG.items(): 
    temp_new_column = {}
    
    for column in cell_df.columns:
        temp_column = column.split("_")
        temp_header = temp_column[4] + '_' + temp_column[5] + '_result'
        
        temp_new_column[column] = temp_header  # Store the mapping of old name to new name
    
    cell_df.rename(columns=temp_new_column, inplace=True)  # Rename the columns using the dictionary

print ('04-3 Rename header of df in dict')

# 04-4 sum the value in DEG
for cell, cell_df in dict_DEG.items():
    # calculate sum of each column
    temp_sum_value = cell_df.sum(axis=0)
    # add a new row for sum
    """
    The unit of result is Log2(Fold Chang) per gene
    For example, the value of 1.5 indicate that 2^1.5 fold change of expression per gene (average) remained after drug treatment. 
    """
    cell_df.loc['sum'] = temp_sum_value/gene_len

print ('04-4 Calculate sum of each compound')
        
#%% 05
# 05 output the mean value of DEG_sum value for each compound
output = pd.DataFrame()
for cell, cell_df in dict_DEG.items():
    temp_output = pd.DataFrame()
    
    temp_sum_value = pd.DataFrame(cell_df.loc['sum'])
    temp_sum_value.rename(columns = {'sum':cell}, inplace=True)
    temp_output = pd.concat([temp_output, temp_sum_value])
    output = pd.concat([output, temp_output])
print ('Preparing output file')

# merge duplicate index
output = output.groupby(output.index).sum()

# replace 0 with NaN in the entire DataFrame
output.replace(0, np.nan, inplace=True)

# calcualte mean value for each row
output['mean'] = output.mean(axis=1)
output = output.sort_values(by='mean')
print ('Processing output file')

# Get the current date and time
date_current = datetime.datetime.now()
date = date_current.strftime('%Y%m%d') # format the date as a string yyyymmdd

# export
file_name = '02_DEG_output_' + f'{date}' + '.csv'
output.to_csv(file_name)
print ('05 Complete output')
