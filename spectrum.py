
import os
import re
import pandas as pd

current_dir = os.getcwd()
parent_dir = os.path.dirname(current_dir)
files_spectrum='data/input_smiles/spectrum/MassBank-data-2024.06'
root_folder=os.path.join(parent_dir, files_spectrum)

# 正则表达式模式
smiles_pattern = re.compile(r'^CH\$SMILES:\s*(.*)')
annotation_pattern = re.compile(r'^PK\$ANNOTATION:\s*m/z\s+tentative_formula\s+mass_error\(ppm\)')
peak_pattern = re.compile(r'^\s*(\S+)\s+(\S+)\s+(\S+)')

# 初始化DataFrame列表
df_list = []

# 遍历文件夹及其子文件夹中的所有txt文件
for root, dirs, files in os.walk(root_folder):
    for file in files:
        if file.endswith('.txt'):
            file_path = os.path.join(root, file)

            # 初始化变量
            smiles = None
            annotations = []

            # 读取txt文件
            with open(file_path, 'r') as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    # 提取CH$SMILES
                    smiles_match = smiles_pattern.match(line)
                    if smiles_match:
                        smiles = smiles_match.group(1)

                    # 提取PK$ANNOTATION
                    if annotation_pattern.match(line):
                        j = i + 1
                        while j < len(lines) and peak_pattern.match(lines[j]):
                            peak_match = peak_pattern.match(lines[j])
                            if peak_match:
                                mz, formula, error = peak_match.groups()
                                annotations.append({
                                    'm/z': mz,
                                    'tentative_formula': formula,
                                    'mass_error(ppm)': error
                                })
                            j += 1

            # 将结果存储到DataFrame
            if annotations:
                df = pd.DataFrame(annotations)
                df.insert(0, 'CH$SMILES', smiles)
                df['File'] = file  # 添加文件名
                df['Path'] = file_path  # 添加文件路径
                df_list.append(df)

# 合并所有DataFrame
final_df = pd.concat(df_list, ignore_index=True)

# 保存到Excel文件
#output_file = '/mnt/work/code/tian/smiles/data/result/spectrum/CH_SMILES_and_PK_ANNOTATION_All.xlsx'
files_spectrum_out='data/input_smiles/spectrum/CH_SMILES_and_PK_ANNOTATION_All.xlsx'
output_file=os.path.join(parent_dir, files_spectrum_out)
final_df.to_excel(output_file, index=False)

print(f"数据已成功保存到 {output_file}")