import pandas as pd
import os
from find_functional_group import identify_functional_groups
from word_img import produce_word_image

current_dir = os.getcwd()
parent_dir = os.path.dirname(current_dir)
data_files_smiles='data/result/statistics_functional_1328'
foldername=os.path.join(parent_dir, data_files_smiles)

#del o row
# df = pd.read_excel((os.path.join(foldername, 'frag_smiles_main.xlsx')),header=0,index_col=0)
# df = df.loc[~(df == 0).all(axis=1)]
# df = df.dropna(how='all')
# df=df.reset_index(drop=False)
# df.to_excel((os.path.join(foldername, 'frag_smiles_main_deal.xlsx')), index=True)



# # find functional group
# file_path = os.path.join(foldername, 'frag_smiles_main_deal.xlsx') # 替换为你的文件路径
# df = pd.read_excel(file_path)
# #df.set_index('smiles',inplace=True)
#
# for i in range(len(df)):
#     smiles=df.at[i,'smiles']
#     frag=df.at[i,'frag']
#     detected_groups = identify_functional_groups(frag,smiles)
#
#     true_keys = [k for item in detected_groups for k, v in item.items() if v is True]
#     if not true_keys:
#         df.at[i, 'functional_grounp'] = df.at[i, 'frag']
#     else:
#         ture_keys_str =','.join(true_keys)
#         df.at[i, 'functional_grounp'] = ture_keys_str
#
# output_file =os.path.join(foldername, 'frag_smiles_main_deal_functional.xlsx')
# df.to_excel(output_file, index=False)


# #devide multiple functional to single functional group
# import pandas as pd
# df = pd.read_excel(os.path.join(foldername, 'frag_smiles_main_deal_functional.xlsx')) # 替换为你的 Excel 文件名和表单名
#
# # 假设需要处理的列名为 'CombinedElements'，并且每个元素是通过逗号分隔的
# split_df = df['functional_grounp'].str.split(',', expand=True).stack()
# # 重置索引并设置列名
# split_df = split_df.reset_index(level=1, drop=True).reset_index(name='fragments')
# # 将结果保存为新的 Excel 文件
# split_df.to_excel(os.path.join(foldername, 'frag_smiles_main_deal_functional_split.xlsx'),index=False)
# #

# #find frequeny
# import pandas as pd
# # 读取 Excel 文件
# file_path = os.path.join(foldername, 'frag_smiles_main_deal_functional_split.xlsx') # 替换为你的文件路径
# df = pd.read_excel(file_path)
# # 假设你要统计的列名是 'ColumnName'
# column_name = 'fragments'  # 替换为你要统计的列名
# frequency = df[column_name].value_counts()
# # 将结果转换为 DataFrame
# frequency_df = frequency.reset_index()
# frequency_df.columns = [column_name, 'Frequency']
# # 将结果保存到新的 Excel 文件中
# output_file =os.path.join(foldername, 'frag_smiles_main_deal_functional_split_frequency.xlsx')
# frequency_df.to_excel(output_file, index=False)

#produce frequency image
file_path = os.path.join(foldername, 'frag_smiles_main_deal_functional_split_frequency.xlsx') # 替换为你的文件路径
df = pd.read_excel(file_path)
produce_word_image(df)

