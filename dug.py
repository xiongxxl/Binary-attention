from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Draw
from rdkit import Chem
# from rdkit.Chem import FragmentOnBonds
# import os
# # 获取当前工作目录
# ratio=0.9
# current_dir = os.getcwd()
# parent_dir = os.path.dirname(current_dir)
# mole_name='c1cf'
# mol_name_suf = f"{mole_name}_{ratio}"
#

import matplotlib.pyplot as plt
import numpy as np
# 生成示例数据
data = np.random.randn(1000)  # 正态分布随机数据
# 创建直方图
plt.hist(data, bins=30, edgecolor='black')
# 添加标题和标签
plt.title('Frequency Histogram')
plt.xlabel('Data')
plt.ylabel('Frequency')
# 显示图形
plt.show()



# folder_name_mol = os.path.join(parent_dir, mol_name_suf)
# print(folder_name_mol)
# import pandas as pd
# # 假设有多个字典
import pandas as pd
# 示例数据
# def is_subset(A, B):
#     return B.issubset(A)
# # 示例
# A = {1, 2, 3, 4, 5}
# B = {2, 3}
# if is_subset(A, B):
#     print("A 包含 B 的所有元素")
# else:
#     print("A 不包含 B 的所有元素")
#

#
# import pandas as pd
# # 示例 DataFrame
# df = pd.DataFrame({
#     'col1': ['1', '2', '3', '4.5', 'not_a_number']
# })
# # 将 'col1' 列转换为数字，如果无法转换则变为 NaN
# df['col1'] = pd.to_numeric(df['col1'], errors='coerce')
# print(df)
#




# import pandas as pd
# # 示例 DataFrame
# df = pd.DataFrame({
#     'col1': ['1', '2', '3', '4.5']
# })
# # 将 'col1' 列转换为浮点数
# df['col1'] = df['col1'].astype(float)  # 或者 int，如果需要整数
# print(df)









# data = {'Name': ['Alice', 'Bob', 'Charlie'], 'Age': [25, 30, 35]}
# df = pd.DataFrame(data)
# # 加入索引列
# df_with_index = df.reset_index()
# print(df_with_index)





# import pandas as pd
#
# # 示例数组
# arr = [10, 15, 18, 22]
# # 创建一个两列的DataFrame
# data = {'A': [[10, 15, 18, 21], [10, 16, 18, 22], [9, 15, 18, 22], [10, 15, 18, 22]],
#         'B': [100, 101, 102, 103]}  # A列是列表，B列是序号
# df = pd.DataFrame(data)
# # 创建新列用于存储标志符，默认为0
# df['Flag'] = 0
#
#
# # 定义比较函数，检查是否完全相同或仅一个元素不同
# def is_similar(arr1, arr2):
#     # 如果长度不同，直接返回False
#     if len(arr1) != len(arr2):
#         return False
#
#     # 计算两数组中不同元素的数量
#     diff_count = sum([1 for x, y in zip(arr1, arr2) if x != y])
#
#     # 如果不同元素的数量为0或1，返回True
#     return diff_count <= 1
#
# # 遍历DataFrame的每一行，与数组进行比较
# for index, row in df.iterrows():
#     # 如果A列中的列表与数组相同或仅一个元素不同
#     if is_similar(row['A'], arr):
#         df.at[index, 'Flag'] = 1  # 将标志符置为1
# # 将满足条件的行及其序号写入Excel文件
# df[df['Flag'] == 1].to_excel('result.xlsx', index=True)
# print("完成！结果已写入 result.xlsx 文件。")






# import pandas as pd
# # 创建示例数组和DataFrame
# arr = [10, 15, 18, 22]
# # 创建一个两列的DataFrame
# data = {'A': [9, 15, 20, 24], 'B': [100, 101, 102, 103]}  # A列和B列
# df = pd.DataFrame(data)
# # 创建新列用于存储标志符，默认为0
# df['Flag'] = 0
# # 定义一个函数用于比较数组中的每个元素和DataFrame的A列
# def compare_and_flag(array, df):
#     # 遍历数组中的每个元素
#     for value in array:
#         # 遍历DataFrame中的每一行，比较A列
#         for index, row in df.iterrows():
#             # 如果A列的值与数组元素相同或相差1
#             if abs(row['A'] - value) <= 1:
#                 df.at[index, 'Flag'] = 1  # 将标志符置为1
#     return df
# # 调用比较函数
# df = compare_and_flag(arr, df)
# # 将满足条件的行写入Excel文件，包括序号
# df[df['Flag'] == 1].to_excel('result.xlsx', index=True)
# print("完成！结果已写入 result.xlsx 文件。")




# data1 = {
#     "Column1": [1, 2, 3],
#     "Column2": ['A', 'B', 'C']
# }
# data2 = {
#     "Column1": [4, 5, 6],
#     "Column2": ['D', 'E', 'F']
# }
# data3 = {
#     "Column1": [7, 8, 9],
#     "Column2": ['G', 'H', 'I']
# }
# # 创建一个空的 DataFrame 来存储结果
# df_combined = pd.DataFrame()
# # 使用 pd.concat() 代替 append()
# df_combined = pd.concat([df_combined, pd.DataFrame(data1)], ignore_index=True)
# df_combined = pd.concat([df_combined, pd.DataFrame(data2)], ignore_index=True)
# df_combined = pd.concat([df_combined, pd.DataFrame(data3)], ignore_index=True)
# # 将结果保存到 Excel 文件
# df_combined.to_excel("combined_dicts_concat_output.xlsx", index=False)



# import pandas as pd
# # 假设有多个字典
# data1 = {
#     "Column1": [1, 2, 3],
#     "Column2": ['A', 'B', 'C']
# }
# data2 = {
#     "Column1": [4, 5, 6],
#     "Column2": ['D', 'E', 'F']
# }
# data3 = {
#     "Column1": [7, 8, 9],
#     "Column2": ['G', 'H', 'I']
# }
# # 创建一个空的 DataFrame 来存储结果
# df_combined = pd.DataFrame()
# # 依次将字典转换为 DataFrame 并使用 append() 添加
# df_combined = df_combined.append(pd.DataFrame(data1), ignore_index=True)
# df_combined = df_combined.append(pd.DataFrame(data2), ignore_index=True)
# df_combined = df_combined.append(pd.DataFrame(data3), ignore_index=True)
# # 将结果保存到 Excel 文件
# df_combined.to_excel("combined_dicts_append_output.xlsx", index=False)





# import pandas as pd
# from wordcloud import WordCloud
# import matplotlib.pyplot as plt
# # 读取 Excel 文件中的某一列数据
# df = pd.read_excel('frag_smiles_main_deal_frequency.xlsx', sheet_name='Sheet1')  # 替换为你的 Excel 文件名和表单名
# text_column = df['frag']  # 替换为你的目标列名
# # 将列中的所有文本连接成一个大的字符串
# text = " ".join(text_column.astype(str).tolist())
# # 生成词云
# wordcloud = WordCloud(width=800, height=400, background_color='white').generate(text)
# # 显示词云
# plt.figure(figsize=(10, 5))
# plt.imshow(wordcloud, interpolation="bilinear")
# plt.axis("off")
# plt.show()
# # 保存词云图像（可选）
# wordcloud.to_file("excel_column_wordcloud.png")
#
#
#
#
# from wordcloud import WordCloud
# import matplotlib.pyplot as plt
# # 输入的一句话
# sentence = "The quick brown fox jumps over the lazy dog."
# # 生成词云
# wordcloud = WordCloud(width=800, height=400, background_color='white').generate(sentence)
# # 显示词云
# plt.figure(figsize=(10, 5))
# plt.imshow(wordcloud, interpolation="bilinear")
# plt.axis("off")
# plt.show()
#
#
#
# # 打开TXT文件以写入
# with open('multiple_dicts.txt', 'w') as file:
#     for i, d in enumerate(dict_list):
#         # 写入字典标题
#         file.write(f'Dictionary {i+1}:\n')
#         # 写入字典内容
#         for key, value in d.items():
#             file.write(f'  {key}: {value}\n')
#         # 添加换行符以分隔字典
#         file.write('\n')
#
# print("字典已成功写入TXT文件")


# current_dir = os.getcwd()
# print(f"当前目录: {current_dir}")
# # 倒退到上一个文件夹地址
# parent_dir = os.path.dirname(current_dir)
# print(f"上一级目录: {parent_dir}")
# # 进入上一个文件夹
# os.chdir(parent_dir)
# new_current_dir = os.getcwd()
# print(f"新的当前目录: {new_current_dir}")











# from openpyxl import load_workbook
# # 指定已有的Excel文件路径
# excel_path = 'example.xlsx'
# # 加载已有的Excel工作簿
# wb = load_workbook(filename=excel_path)
# # 选择要写入数据的工作表，这里假设是第一个工作表
# sheet = wb.active
# # 假设有一个列表，包含要写入的多行数据
# data_to_append = [
#     ["姓名", "年龄", "性别"],
#     ["赵六", 25, "男"],
#     ["孙七", 35, "女"]
# ]
# # 确定写入数据的起始行号，这里假设在现有数据之后写入
# start_row = sheet.max_row + 1
# # 逐行写入数据
# for row_data in data_to_append:
#     sheet.append(row_data)
# # 保存Excel文件，覆盖原有文件
# wb.save(filename=excel_path)
#
#
#
#
#
# from openpyxl import Workbook, load_workbook
# # 示例数据列表
# data = [
#     ["Name", "Age", "City"],
#     ["Alice", 30, "New York"],
#     ["Bob", 25, "Los Angeles"],
#     ["Charlie", 35, "Chicago"]
# ]
# # 创建一个新的Excel工作簿和工作表，或者加载现有的Excel文件
# file_path = 'example.xlsx'
# try:
#     # 尝试加载现有的Excel文件
#     wb = load_workbook(file_path)
#     ws = wb.active
# except FileNotFoundError:
#     # 如果文件不存在，则创建一个新的工作簿和工作表
#     wb = Workbook()
#     ws = wb.active
#     ws.title = "Sheet1"
# # 将数据逐行写入工作表
# for row in data:
#     ws.append(row)
# # 保存Excel文件
# wb.save(file_path)
# print(f'已将数据保存到{file_path}')



# from rdkit import Chem
# from rdkit.Chem import Descriptors, rdMolDescriptors
# import pandas as pd
# # 定义你想检测的官能团SMARTS模式
# functional_groups = {
#     'Hydroxyl': '[OX2H]',
#     'Amine': '[NX3;H2,H1;!$(NC=O)]',
#     'Carboxyl': 'C(=O)[OH]',
#     'Aldehyde': '[CX3H1](=O)[#6]',
#     'Ketone': '[CX3](=O)[#6]',
#     'Ester': '[CX3](=O)[OX2H0][#6]',
#     'Ether': '[OD2]([#6])[#6]',
# }
# # 示例分子列表
# smiles_list = [
#     'CCO',  # 乙醇
#     'CC(=O)O',  # 醋酸
#     'CC(=O)C',  # 丙酮
#     'CC(=O)OC',  # 甲酸乙酯
# ]
# # 创建一个空的DataFrame
# columns = ['SMILES'] + list(functional_groups.keys())
# df = pd.DataFrame(columns=columns)
# # 遍历每个分子并检测官能团
# for smiles in smiles_list:
#     mol = Chem.MolFromSmiles(smiles)
#     row = {'SMILES': smiles}
#     for group_name, smarts in functional_groups.items():
#         patt = Chem.MolFromSmarts(smarts)
#         matches = mol.GetSubstructMatches(patt)
#         row[group_name] = len(matches)
#     df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
# # 将DataFrame保存为Excel文件d
# df.to_excel('functional_groups.xlsx', index=False)

#
# from rdkit import Chem
# from rdkit.Chem import rdMolDescriptors
# from rdkit.Chem import Draw
# import pandas as pd
# # 定义要检测的官能团（示例中的一些常见官能团的 SMARTS 模式）
# functional_groups = {
#     'hydroxyl': '[OX2H]',  # 羟基
#     'carbonyl': '[CX3]=[OX1]',  # 羰基
#     'amine': '[NX3;H2,H1;!$(NC=O)]',  # 胺基
#     'carboxyl': '[CX3](=O)[OX2H1]',  # 羧基
# }
# # 读取分子（这里使用一些示例的 SMILES）
# smiles_list = ["CCO", "CC(=O)O", "CCN", "CC(C(=O)O)O"]
# molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]


# # 创建数据框来保存结果
# data = []
# # 检测每个分子中的官能团
# for mol in molecules:
#     mol_data = {'smiles': Chem.MolToSmiles(mol)}
#     for group_name, smarts in functional_groups.items():
#         pattern = Chem.MolFromSmarts(smarts)
#         mol_data[group_name] = mol.HasSubstructMatch(pattern)
#     data.append(mol_data)
# # 创建 pandas 数据框
# df = pd.DataFrame(data)
# # 保存到 Excel 文件
# df.to_excel("functional_groups.xlsx", index=False)
# print("官能团检测结果已保存到 functional_groups.xlsx")


# def split_string_by_dot(input_string):
#     # 使用split方法按照'.'分割字符串
#     elements = input_string.split('.')
#     return elements
#
# # 示例输入
# input_string = "H2O.CC.O2"
# elements = split_string_by_dot(input_string)
# print("Split Elements:", elements)

# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     if mol is None:
#         raise ValueError("Invalid SMILES string")
#
#     # 标记需要保留的原子
#     atoms_to_keep = set(atom_indices)
#
#     # 创建一个新的分子用于存放提取的部分
#     new_mol = Chem.RWMol()
#     atom_map = {}
#
#     for atom in mol.GetAtoms():
#         if atom.GetIdx() in atoms_to_keep:
#             new_atom_idx = new_mol.AddAtom(atom)
#             atom_map[atom.GetIdx()] = new_atom_idx
#
#     # 复制保留原子之间的键
#     for bond in mol.GetBonds():
#         begin_idx = bond.GetBeginAtomIdx()
#         end_idx = bond.GetEndAtomIdx()
#         if begin_idx in atoms_to_keep and end_idx in atoms_to_keep:
#             new_mol.AddBond(atom_map[begin_idx], atom_map[end_idx], bond.GetBondType())
#
#     # 将新分子转为Mol对象，去芳香化避免错误
#     new_mol = new_mol.GetMol()
#     Chem.Kekulize(new_mol, clearAromaticFlags=True)
#
#     # 获取独立的片段
#     frags = rdmolops.GetMolFrags(new_mol, asMols=True)
#
#     # 获取每个片段的SMILES表示
#     frag_smiles_list = [Chem.MolToSmiles(frag) for frag in frags]
#
#     return frag_smiles_list
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 2, 7, 8, 9]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragments = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragments SMILES:", extracted_fragments)

# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     if mol is None:
#         raise ValueError("Invalid SMILES string")
#
#     # 标记需要保留的原子
#     atoms_to_keep = set(atom_indices)
#
#     # 创建一个新的分子用于存放提取的部分
#     new_mol = Chem.RWMol()
#     atom_map = {}
#
#     for atom in mol.GetAtoms():
#         if atom.GetIdx() in atoms_to_keep:
#             new_atom_idx = new_mol.AddAtom(atom)
#             atom_map[atom.GetIdx()] = new_atom_idx
#
#     # 复制保留原子之间的键
#     for bond in mol.GetBonds():
#         begin_idx = bond.GetBeginAtomIdx()
#         end_idx = bond.GetEndAtomIdx()
#         if begin_idx in atoms_to_keep and end_idx in atoms_to_keep:
#             new_mol.AddBond(atom_map[begin_idx], atom_map[end_idx], bond.GetBondType())
#
#     # 将新分子转为Mol对象，避免芳香化问题
#     new_mol = new_mol.GetMol()
#     Chem.SanitizeMol(new_mol)
#
#     # 获取独立的片段
#     frags = rdmolops.GetMolFrags(new_mol, asMols=True)
#
#     # 获取每个片段的SMILES表示
#     frag_smiles_list = [Chem.MolToSmiles(frag) for frag in frags]
#
#     return frag_smiles_list
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 2, 7, 8, 9]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragments = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragments SMILES:", extracted_fragments)

# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 创建一个新的分子编辑对象
#     rw_mol = Chem.RWMol(mol)
#
#     # 标记需要保留的原子
#     atoms_to_keep = set(atom_indices)
#
#     # 创建一个新的分子用于存放提取的部分
#     new_mol = Chem.RWMol()
#     atom_map = {}
#
#     for atom in rw_mol.GetAtoms():
#         if atom.GetIdx() in atoms_to_keep:
#             new_atom_idx = new_mol.AddAtom(atom)
#             atom_map[atom.GetIdx()] = new_atom_idx
#
#     # 复制保留原子之间的键
#     for bond in rw_mol.GetBonds():
#         begin_idx = bond.GetBeginAtomIdx()
#         end_idx = bond.GetEndAtomIdx()
#         if begin_idx in atoms_to_keep and end_idx in atoms_to_keep:
#             new_mol.AddBond(atom_map[begin_idx], atom_map[end_idx], bond.GetBondType())
#
#     # 获取独立的片段
#     frags = rdmolops.GetMolFrags(new_mol, asMols=True)
#
#     # 获取每个片段的SMILES表示
#     frag_smiles_list = [Chem.MolToSmiles(frag) for frag in frags]
#
#     return frag_smiles_list
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 2, 7, 8, 9]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragments = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragments SMILES:", extracted_fragments)

# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 创建一个新的分子编辑对象
#     rw_mol = Chem.RWMol(mol)
#
#     # 标记需要删除的原子
#     atoms_to_remove = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetIdx() not in atom_indices]
#
#     # 逆序删除原子，以防索引改变
#     for atom_idx in sorted(atoms_to_remove, reverse=True):
#         rw_mol.RemoveAtom(atom_idx)
#
#     # 将分子转回普通分子对象
#     frag_mol = rw_mol.GetMol()
#
#     # 获取片段的SMILES表示
#     frag_smiles = Chem.MolToSmiles(frag_mol)
#
#     return frag_smiles
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 2, 7, 8, 9]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragment_smiles = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragment SMILES:", extracted_fragment_smiles)

# def is_single_atom(smiles):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#     # 如果解析失败，则返回False
#     if mol is None:
#         return False
#     # 返回分子中原子的数量是否为1
#     return mol.GetNumAtoms() == 1
#
# # 示例输入
# smiles_list = [
#     "C",      # 单个碳原子
#     "CC",     # 乙烷
#     "O",      # 单个氧原子
#     "H2O",    # 水（无效SMILES）
#     "N",      # 单个氮原子
#     "CCO"     # 乙醇
# ]
#
# # 检查每个SMILES字符串
# for smiles in smiles_list:
#     print(f"SMILES: {smiles} -> Is single atom: {is_single_atom(smiles)}")
















# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 创建一个新的分子编辑对象
#     rw_mol = Chem.RWMol(mol)
#
#     # 标记需要删除的原子
#     atoms_to_remove = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetIdx() not in atom_indices]
#
#     # 逆序删除原子，以防索引改变
#     for atom_idx in sorted(atoms_to_remove, reverse=True):
#         rw_mol.RemoveAtom(atom_idx)
#
#     # 将分子转回普通分子对象
#     frag_mol = rw_mol.GetMol()
#
#     # 获取片段的SMILES表示
#     frag_smiles = Chem.MolToSmiles(frag_mol)
#
#     return frag_smiles

#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0,1,3]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragment_smiles = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragment SMILES:", extracted_fragment_smiles)
#
# #
# def extract_specified_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 创建一个新的分子编辑对象
#     rw_mol = Chem.RWMol(mol)
#
#     # 标记需要保留的原子
#     atoms_to_keep = set(atom_indices)
#
#     # 删除不在原子索引列表中的原子
#     for atom in reversed(rw_mol.GetAtoms()):
#         if atom.GetIdx() not in atoms_to_keep:
#             rw_mol.RemoveAtom(atom.GetIdx())
#
#     # 将分子转回普通分子对象
#     frag_mol = rw_mol.GetMol()
#
#     # 获取片段的SMILES表示
#     frag_smiles = Chem.MolToSmiles(frag_mol)
#
#     return frag_smiles
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 2, 7, 8, 9]  # 指定要保留的原子位置
#
# # 提取片段并输出结果
# extracted_fragment_smiles = extract_specified_atoms(smiles, atom_indices)
# print("Extracted Fragment SMILES:", extracted_fragment_smiles)


# def split_molecule_by_atoms(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 转换为分子编辑对象
#     rw_mol = Chem.RWMol(mol)
#
#     # 切割指定原子位置
#     bonds_to_cut = []
#     for i in atom_indices:
#         atom = rw_mol.GetAtomWithIdx(i)
#         for bond in atom.GetBonds():
#             bonds_to_cut.append(bond.GetIdx())
#
#     # 切割分子
#     fragments = rdmolops.FragmentOnBonds(rw_mol, bonds_to_cut, addDummies=False)
#
#     # 获取片段的SMILES表示
#     frags = Chem.GetMolFrags(fragments, asMols=True)
#     frags_smiles = [Chem.MolToSmiles(frag) for frag in frags]
#
#     return frags_smiles
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [1, 6]  # 指定要切割的原子位置
#
# # 切割分子并输出结果
# fragments_smiles = split_molecule_by_atoms(smiles, atom_indices)
# for i, frag_smiles in enumerate(fragments_smiles):
#     print(f"Fragment {i + 1}: {frag_smiles}")
#
#
# def split_molecule(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 获取要切割的键
#     bonds_to_cut = []
#     for i in atom_indices:
#         atom = mol.GetAtomWithIdx(i)
#         for bond in atom.GetBonds():
#             bond_idx = bond.GetIdx()
#             if bond_idx not in bonds_to_cut:
#                 bonds_to_cut.append(bond_idx)
#
#     # 切割分子
#     fragments = FragmentOnBonds(mol, bonds_to_cut)
#
#     # 获取片段的SMILES表示
#     frags_smiles = Chem.MolToSmiles(fragments)
#
#     return frags_smiles
#
#
# # 示例输入：SMILES字符串和原子位置
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 阿司匹林
# atom_indices = [0, 1, 10]  # 指定要切割的原子位置
#
# # 切割分子并输出结果
# fragments_smiles = split_molecule(smiles, atom_indices)
# print("Fragmented SMILES:", fragments_smiles)

#
# # 定义常见官能团的SMARTS模式及其化学式
# functional_groups = {
#     "Alkane": {"smarts": "[CX4]", "formula": "C-H"},
#     "Alkene": {"smarts": "[CX3]=[CX3]", "formula": "C=C"},
#     "Alkyne": {"smarts": "[CX2]#[CX2]", "formula": "C≡C"},
#     "Aromatic": {"smarts": "c1ccccc1", "formula": "C6H5"},
#     "Halide": {"smarts": "[F,Cl,Br,I]", "formula": "X (F, Cl, Br, I)"},
#     "Alcohol": {"smarts": "[OX2H]", "formula": "R-OH"},
#     "Phenol": {"smarts": "c1ccc(O)cc1", "formula": "C6H5OH"},
#     "Ether": {"smarts": "[OD2]([#6])[#6]", "formula": "R-O-R"},
#     "Aldehyde": {"smarts": "[CX3H1](=O)[#6]", "formula": "R-CHO"},
#     "Ketone": {"smarts": "[CX3](=O)[#6]", "formula": "R-CO-R"},
#     "Carboxylic Acid": {"smarts": "[CX3](=O)[OX2H1]", "formula": "R-COOH"},
#     "Ester": {"smarts": "[CX3](=O)[OX2][#6]", "formula": "R-COO-R"},
#     "Amide": {"smarts": "[NX3][CX3](=[OX1])[#6]", "formula": "R-CONH2"},
#     "Amine": {"smarts": "[NX3][#6]", "formula": "R-NH2"},
#     "Nitrate": {"smarts": "[NX3](=O)([OX1-])[OX1-]", "formula": "R-NO3"},
#     "Nitro": {"smarts": "[NX3](=O)[OX1-]", "formula": "R-NO2"},
#     "Sulfonic Acid": {"smarts": "S(=O)(=O)[O-]", "formula": "R-SO3H"},
#     "Thiol": {"smarts": "[SX2H]", "formula": "R-SH"},
#     "Thioether": {"smarts": "[SX2][#6]", "formula": "R-S-R"},
#     "Disulfide": {"smarts": "[SX2][SX2]", "formula": "R-S-S-R"},
#     "Sulfoxide": {"smarts": "[SX3](=O)[#6]", "formula": "R-S(=O)-R"},
#     "Sulfone": {"smarts": "[SX4](=O)(=O)[#6]", "formula": "R-SO2-R"},
#     "Phosphine": {"smarts": "[PX3]", "formula": "R3P"},
#     "Phosphate": {"smarts": "P(=O)(O)(O)O", "formula": "R-O-PO3H2"},
#     "Isocyanate": {"smarts": "N=C=O", "formula": "R-N=C=O"},
#     "Isothiocyanate": {"smarts": "N=C=S", "formula": "R-N=C=S"}
# }
#
#
# def identify_functional_groups(smiles):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 存储检测到的官能团及其化学式
#     detected_groups = []
#
#     # 查找官能团
#     for group_name, properties in functional_groups.items():
#         pattern = Chem.MolFromSmarts(properties["smarts"])
#         if mol.HasSubstructMatch(pattern):
#             detected_groups.append((group_name, properties["formula"]))
#
#     return detected_groups
#
#
# # 示例输入：SMILES字符串
# smiles = "CC(=O)O"  # 乙酸
#
# # 识别并输出官能团
# detected_groups = identify_functional_groups(smiles)
# for group_name, formula in detected_groups.items():
#     print(f"Detected functional group: {group_name}, Chemical formula: {formula}")

# # 定义常见官能团的SMARTS模式
# functional_groups = {
#     "Alkane": "[CX4]",
#     "Alkene": "[CX3]=[CX3]",
#     "Alkyne": "[CX2]#[CX2]",
#     "Aromatic": "c1ccccc1",
#     "Halide": "[F,Cl,Br,I]",
#     "Alcohol": "[OX2H]",
#     "Phenol": "c1ccc(O)cc1",
#     "Ether": "[OD2]([#6])[#6]",
#     "Aldehyde": "[CX3H1](=O)[#6]",
#     "Ketone": "[CX3](=O)[#6]",
#     "Carboxylic Acid": "[CX3](=O)[OX2H1]",
#     "Ester": "[CX3](=O)[OX2][#6]",
#     "Amide": "[NX3][CX3](=[OX1])[#6]",
#     "Amine": "[NX3][#6]",
#     "Nitrate": "[NX3](=O)([OX1-])[OX1-]",
#     "Nitro": "[NX3](=O)[OX1-]",
#     "Sulfonic Acid": "S(=O)(=O)[O-]",
#     "Thiol": "[SX2H]",
#     "Thioether": "[SX2][#6]",
#     "Disulfide": "[SX2][SX2]",
#     "Sulfoxide": "[SX3](=O)[#6]",
#     "Sulfone": "[SX4](=O)(=O)[#6]",
#     "Phosphine": "[PX3]",
#     "Phosphate": "P(=O)(O)(O)O",
#     "Isocyanate": "N=C=O",
#     "Isothiocyanate": "N=C=S"
# }
#
#
# def identify_functional_groups(smiles):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 存储检测到的官能团
#     detected_groups = set()
#
#     # 查找官能团
#     for group_name, smarts in functional_groups.items():
#         pattern = Chem.MolFromSmarts(smarts)
#         if mol.HasSubstructMatch(pattern):
#             detected_groups.add(group_name)
#
#     return detected_groups
#
#
# # 示例输入：SMILES字符串
# #smiles = "CC(=O)O"  # 乙酸
# smiles="NC(=S)Nc1ccccc1"
# # 识别并输出官能团
# detected_groups = identify_functional_groups(smiles)
# print("Detected functional groups:", detected_groups)

## function group chinese
# Alkane: "[CX4]" 匹配所有四价的碳原子，即饱和碳（烷烃）。
# Alkene: "[CX3]=[CX3]" 匹配双键碳原子（烯烃）。
# Alkyne: "[CX2]#[CX2]" 匹配三键碳原子（炔烃）。
# Aromatic: "c1ccccc1" 匹配苯环。
# Halide: "[F,Cl,Br,I]" 匹配卤素原子（氟、氯、溴、碘）。
# Alcohol: "[OX2H]" 匹配羟基（酒精）。
# Phenol: "c1ccc(O)cc1" 匹配苯酚。
# Ether: "OD2[#6]" 匹配醚。
# Aldehyde: "CX3H1[#6]" 匹配醛。
# Ketone: "CX3[#6]" 匹配酮。
# Carboxylic Acid: "CX3[OX2H1]" 匹配羧酸。
# Ester: "CX3[OX2][#6]" 匹配酯。
# Amide: "[NX3]CX3[#6]" 匹配酰胺。
# Amine: "[NX3][#6]" 匹配胺。
# Nitrate: "NX3([OX1-])[OX1-]" 匹配硝酸盐。
# Nitro: "NX3[OX1-]" 匹配硝基。
# Sulfonic Acid: "S(=O)(=O)[O-]" 匹配磺酸。
# Thiol: "[SX2H]" 匹配硫醇。
# Thioether: "[SX2][#6]" 匹配硫醚。
# Disulfide: "[SX2][SX2]" 匹配二硫化物。
# Sulfoxide: "SX3[#6]" 匹配亚砜。
# Sulfone: "SX4(=O)[#6]" 匹配砜。
# Phosphine: "[PX3]" 匹配膦。
# Phosphate: "P(=O)(O)(O)O" 匹配磷酸。
# Isocyanate: "N=C=O" 匹配异氰酸酯。
# Isothiocyanate: "N=C=S" 匹配异硫氰酸酯。




#
# # 定义常见官能团的SMARTS模式
# functional_groups = {
#     "Hydroxyl": "[O]",
#     "Carbonyl": "[C=O]",
#     "Carboxyl": "[C(=O)O]",
#     "Amino": "[NH2]",
#     "Phenyl": "c1ccccc1",
#     "Ether": "[C-O-C]",
#     "Ester": "[C(=O)O-C]",
#     "Amide": "[C(=O)N]",
#     "Nitrate": "[N+](=O)[O-]"
# }
#
#
# def identify_functional_groups(smiles):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 存储检测到的官能团
#     detected_groups = set()
#
#     # 查找官能团
#     for group_name, smarts in functional_groups.items():
#         pattern = Chem.MolFromSmarts(smarts)
#         if mol.HasSubstructMatch(pattern):
#             detected_groups.add(group_name)
#
#     return detected_groups
#
#
# # 示例输入：SMILES字符串
# smiles = "CC(=O)O"  # 乙酸
#
# # 识别并输出官能团
# detected_groups = identify_functional_groups(smiles)
# print("Detected functional groups:", detected_groups)
#
# # from rdkit import Chem
# #
# # # 定义常见官能团的SMARTS模式
# # functional_groups = {
# #     "Hydroxyl": "[OH]",
# #     "Carbonyl": "[C=O]",
# #     "Carboxyl": "[C(=O)O]",
# #     "Amino": "[NH2]",
# #     "Phenyl": "c1ccccc1",
# #     "Ether": "[C-O-C]",
# #     "Ester": "[C(=O)O-C]",
# #     "Amide": "[C(=O)N]",
# #     "Nitrate": "[N+](=O)[O-]"
# # }
# #
# #
# # def check_functional_groups(smiles, atom_indices):
# #     # 使用RDKit解析SMILES字符串
# #     mol = Chem.MolFromSmiles(smiles)
# #
#     # 存储每个原子是否属于某个官能团
#     atom_functional_groups = {idx: [] for idx in atom_indices}
#
#     # 查找官能团
#     for group_name, smarts in functional_groups.items():
#         pattern = Chem.MolFromSmarts(smarts)
#         matches = mol.GetSubstructMatches(pattern)
#         for match in matches:
#             for idx in match:
#                 if idx in atom_indices:
#                     atom_functional_groups[idx].append(group_name)
#
#     return atom_functional_groups
#
#
# # 示例输入：SMILES字符串和对应的原子序号列表
# smiles = "CC(=O)O"  # 乙酸
# atom_indices = [0, 1]  # 原子序号列表
#
# # 检查这些原子是否属于某个官能团
# result = check_functional_groups(smiles, atom_indices)
#
# # 输出结果
# for atom_idx, groups in result.items():
#     print(f"Atom {atom_idx}: {groups}")




# from rdkit import Chem
#
# # 定义官能团SMARTS模式
# functional_groups = {
#     "Hydroxyl": "[OH]",
#     "Carbonyl": "[C=O]",
#     "Carboxyl": "[C(=O)O]",
#     "Amino": "[NH2]",
#     "Phenyl": "c1ccccc1",
#     "Ether": "[C-O-C]",
#     "Ester": "[C(=O)O-C]",
#     "Amide": "[C(=O)N]",
#     "Nitrate": "[N+](=O)[O-]"
# }
#
#
# def check_functional_groups(smiles, atom_indices):
#     # 使用RDKit解析SMILES字符串
#     mol = Chem.MolFromSmiles(smiles)
#
#     # 存储每个原子是否属于某个官能团
#     atom_functional_groups = {idx: [] for idx in atom_indices}
#
#     # 查找官能团
#     for group_name, smarts in functional_groups.items():
#         pattern = Chem.MolFromSmarts(smarts)
#         matches = mol.GetSubstructMatches(pattern)
#         for match in matches:
#             for idx in match:
#                 if idx in atom_indices:
#                     atom_functional_groups[idx].append(group_name)
#
#     return atom_functional_groups
#
#
# # 输入SMILES字符串和对应的原子序号列表
# smiles = "CC(=O)O"
# atom_indices = [1, 2, 3, 4]  # 例子中的原子序号
#
# # 检查这些原子是否属于某个官能团
# result = check_functional_groups(smiles, atom_indices)
#
# # 输出结果
# for atom_idx, groups in result.items():
#     print(f"Atom {atom_idx}: {groups}")

#
# from rdkit import Chem
# from rdkit.Chem import rdMolDescriptors
#
# # 输入SMILES字符串和对应的原子序号列表
# smiles = "CCO"
# atom_indices = [1, 2, 3]  # 举例的原子序号
#
# # 使用RDKit解析SMILES字符串
# mol = Chem.MolFromSmiles(smiles)
#
# # 定义常见官能团SMARTS模式
# functional_groups = {
#     "Hydroxyl": "[OH]",
#     "Carbonyl": "[C=O]",
#     "Carboxyl": "[C(=O)O]",
#     "Amino": "[NH2]",
#     "Phenyl": "c1ccccc1",
#     "Ether": "[C-O-C]",
#     "Ester": "[C(=O)O-C]",
#     "Amide": "[C(=O)N]",
#     "Nitrate": "[N+](=O)[O-]"
# }
#
# # 查找官能团
# detected_groups = []
# for group_name, smarts in functional_groups.items():
#     pattern = Chem.MolFromSmarts(smarts)
#     if mol.HasSubstructMatch(pattern):
#         detected_groups.append(group_name)
#
# # 输出存在的官能团
# print("Detected functional groups:", detected_groups)

#
# from rdkit import Chem
# from rdkit.Chem import Draw
#
# def find_functional_groups(smiles, functional_groups):
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         return None
#
#     results = {}
#     for name, smarts in functional_groups.items():
#         patt = Chem.MolFromSmarts(smarts)
#         matches = mol.GetSubstructMatches(patt)
#         if matches:
#             results[name] = matches
#
#     return results
#
# # 定义官能团的SMARTS
# functional_groups = {
#     "Hydroxyl": "[OX2H]",        # 羟基
#     "Carbonyl": "[CX3]=[OX1]",   # 羰基
#     "Amine": "[NX3;H2,H1;!$(NC=O)]", # 胺基
#     "Carboxyl": "C(=O)[OH]",     # 羧基
#     # 可以添加更多官能团
# }
#
# # 测试示例
# smiles = "CC(=O)Oc1ccccc1C(=O)O"  # 乙酸苯酯
#
# results = find_functional_groups(smiles, functional_groups)
#
# if results:
#     for fg, matches in results.items():
#         print(f"{fg} found at positions: {matches}")
# else:
#     print("No functional groups found.")
#
# # 绘制分子结构并标注官能团位置
# mol = Chem.MolFromSmiles(smiles)
# Draw.MolToImage(mol, highlightAtoms=[atom for match in results.values() for group in match for atom in group])
#

# def find_functional_groups(smiles, functional_groups):
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         return None
#
#     results = {}
#     for name, smarts in functional_groups.items():
#         patt = Chem.MolFromSmarts(smarts)
#         matches = mol.GetSubstructMatches(patt)
#         if matches:
#             results[name] = matches
#
#     return results
#
# # 定义官能团的SMARTS
# functional_groups = {
#     "Hydroxyl": "[OX2H]",        # 羟基
#     "Carbonyl": "[CX3]=[OX1]",   # 羰基
#     "Amine": "[NX3;H2,H1;!$(NC=O)]", # 胺基
#     "Carboxyl": "C(=O)[OH]",     # 羧基
#     # 可以添加更多官能团
# }
#
# # 测试示例
# smiles = "CC(=O)Oc1ccccc1C(=O)O"  # 乙酸苯酯
# results = find_functional_groups(smiles, functional_groups)
#
# if results:
#     for fg, matches in results.items():
#         print(f"{fg} found at positions: {matches}")
# else:
#     print("No functional groups found.")



# def find_and_remove_non_letters(input_string):
#     non_letter_positions = []
#     clean_string = ""
#
#     for i, char in enumerate(input_string):
#         if not char.isalpha():
#             non_letter_positions.append(i)
#         else:
#             clean_string += char
#
#     return non_letter_positions, clean_string
#
# # 测试示例
# input_string = "Hello, World! 123"
# positions, clean_string = find_and_remove_non_letters(input_string)
#
# print("非字母字符的位置:", positions)
# print("删除非字母字符后的字符串:", clean_string)


# import re
#
# def find_non_alphanumeric_positions(s):
#     # 使用正则表达式找到所有非字母和数字的元素位置
#     return [(match.start(), match.group()) for match in re.finditer(r'\W', s)]
#
# # 示例字符串
# s = "Hello, World! 123."
#
# positions = find_non_alphanumeric_positions(s)
# print(positions)

# import re
# import numpy as np
#
# def find_non_alphanumeric_positions(s):
#     # 使用正则表达式找到所有非字母和数字的元素位置
#     return [(match.start(), match.group()) for match in re.finditer(r'\W', s)]
#
# def positions_matrix(positions, length):
#     # 创建一个与字符串长度相同的零矩阵
#     matrix = np.zeros((1, length), dtype=int)
#     for pos, char in positions:
#         matrix[0, pos] = pos + 1  # 位置从1开始编号
#     return matrix
#
# # 示例字符串
# s = "Hello, World! 123."
#
# positions = find_non_alphanumeric_positions(s)
# matrix = positions_matrix(positions, len(s))
#
# print("Positions and characters:", positions)
# print("Positions matrix:\n", matrix)

# import re
#
#
# def find_non_alphanumeric_positions_and_remove(s):
#     # 使用正则表达式找到所有非字母和数字的元素位置
#     positions = [match.start() for match in re.finditer(r'\W', s)]
#
#     # 删除非字母和数字元素后生成新的字符串
#     new_string = re.sub(r'\W', '', s)
#
#     return positions, new_string
#
#
# # 示例字符串
# s = "Hello, World! 123."
#
# positions, new_string = find_non_alphanumeric_positions_and_remove(s)
#
# print("Positions of non-alphanumeric characters:", positions)
# print("String after removing non-alphanumeric characters:", new_string)


# import re
#
#
# def find_non_alphanumeric_positions_and_remove(lst):
#     # 找到所有非字母和数字的元素位置
#     positions = [i for i, elem in enumerate(lst) if not re.match(r'\w', str(elem))]
#
#     # 删除非字母和数字元素后生成新的列表
#     new_list = [elem for i, elem in enumerate(lst) if i not in positions]
#
#     return positions, new_list
#
#
# # 示例列表
# lst = ['H', 'e', 'l', 'l', 'o', ',', ' ', 'W', 'o', 'r', 'l', 'd', '!', ' ', '1', '2', '3', '.']
#
# positions, new_list = find_non_alphanumeric_positions_and_remove(lst)
#
# print("Positions of non-alphanumeric characters:", positions)
# print("List after removing non-alphanumeric characters:", new_list)



