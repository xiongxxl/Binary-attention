import numpy as np
import os
from tokenizer import tokenize_smiles
from highlight_num_atoms import highlight_num_atoms
from remove_no_alpha import find_and_remove_non_letters
from remove_no_alpha import find_non_alphanumeric_positions_and_remove
from exact_smile_fragment import extract_specified_atoms
from find_functional_group import identify_functional_groups
import pandas as pd
from PIL import Image
import io

def remove_rows_and_columns(matrix, rows_to_remove, cols_to_remove):
    # 转换输入列表为NumPy数组
    matrix_np = np.array(matrix)
    # 删除指定的行
    matrix_np = np.delete(matrix_np, rows_to_remove, axis=0)
    # 删除指定的列
    matrix_np = np.delete(matrix_np, cols_to_remove, axis=1)
    return matrix_np

#folder_path = './result/attention/first/img/CCCCCCCCCCCCCCO'
# 遍历根目录下的所有文件夹
all_frag_smiles = []
all_frag_functional=[]
def highlight_chemical_atoms (folder_path,filename_without_extension,parent_dir,statistics_fragment_path,statistics_location_path):
    if 'H' in filename_without_extension:
        pass
    else:
        for file_name in os.listdir(folder_path):
            if file_name.endswith('.npy'):
                file_path = os.path.join(folder_path, file_name)
                attn_binary = np.load(file_path)
                #np.fill_diagonal(attn_binary,0) #set diagonal element to 0
                attn_binary_del=attn_binary

                # del non-character
                filename_without_extension_re = tokenize_smiles(filename_without_extension)
                ## delete alpha
                #positions, cleaned_string = find_non_alphanumeric_positions_and_remove(filename_without_extension_re)
                positions, cleaned_string = find_and_remove_non_letters(filename_without_extension_re)
                # np.save(f'{filename_without_extension}.npy',attn_binary_del)
                # print('cleaned_string is:',cleaned_string)
                # print('positions is:', positions)
                # print('attn_binary :', np.shape(attn_binary))
                # print(attn_binary_del)

                # del Delete the corresponding rows and columns
                attn_binary_letter = remove_rows_and_columns(attn_binary_del, positions, positions)
                #print('attn_binary_del:', np.shape(attn_binary_letter))

                # find no_zero_element
                non_zero_positions = np.nonzero(attn_binary_letter)
                non_zero_col=np.unique(non_zero_positions[0]).tolist()
                non_zero_row=np.unique(non_zero_positions[1]).tolist()
                # print("不为零元素的行索引：", non_zero_col)
                # print("不为零元素的列索引：", non_zero_row)

                ## find functonal group
                smiles = filename_without_extension
                atom_indices = non_zero_col
                keep_single_elements,del_single_elements,frag_smiles = extract_specified_atoms(smiles, atom_indices)
                detected_groups = identify_functional_groups(frag_smiles, smiles)

                ## save all_frag_smiles
                file_name_nosuffix = os.path.splitext(file_name)[0]
                frag_smiles_title = {'smiles': smiles}
                frag_smiles_content ={'frag': del_single_elements}
                frag_smiles_aixis = {'axis': file_name_nosuffix}
                all_frag_smiles=frag_smiles_title.copy()
                all_frag_smiles.update(frag_smiles_content)
                all_frag_smiles.update(frag_smiles_aixis)
                ## save to  Excel
                #statistics_excel_path = 'data/result/statistics_functional_single/frag_smiles_main.xlsx'
                folder_frag_excel = os.path.join(parent_dir,statistics_fragment_path)
                df_excel = pd.read_excel(folder_frag_excel)
                df_new = pd.DataFrame(all_frag_smiles)
                df_excel.loc[len(df_excel)] = df_new.iloc[0, :]
                df_excel.to_excel(folder_frag_excel, index=False)

                #save all frag location
                file_name_nosuffix = os.path.splitext(file_name)[0]
                frag_location_title = {'smiles': smiles}
                frag_location_content ={'location': [str(non_zero_col)]}
                frag_location_aixis={'axis': file_name_nosuffix}
                all_frag_location=frag_location_title.copy()
                all_frag_location.update(frag_location_content)
                all_frag_location.update(frag_location_aixis)


               # statistics_excel = 'data/result/statistics_functional_12015/functional_groups_main.xlsx'
                folder_location_excel = os.path.join(parent_dir,statistics_location_path)
                df_excel_location = pd.read_excel(folder_location_excel)
                df_location = pd.DataFrame(all_frag_location)
                df_excel_location.loc[len(df_excel_location)] = df_location.iloc[0, :]
                df_excel_location.to_excel(folder_location_excel, index=False)


                ##highlight atoms
                atom_indices=non_zero_col
                #atom_indices = non_zero_row
                #print(atom_indices)
                img_data = highlight_num_atoms(smiles, atom_indices)
                file_name_no_suffix=os.path.splitext(file_name)[0]
                img_filename=f'{file_name_no_suffix}.jpeg'
                img_path=os.path.join(folder_path,img_filename)
                with open(img_path, "wb") as f:
                    f.write(img_data)
                # img = Image.open(io.BytesIO(img_data))
                # img.show()

        return atom_indices

if __name__ == "__main__":
    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    data_files_npy_simple ='data/result/img_1328/CC(C)(C)c1ccc(-c2ccc([N+](=O)[O-])cc2)cc1'
    folder_path =os.path.join(parent_dir, data_files_npy_simple)
    filename_without_extension='CC(C)(C)c1ccc(-c2ccc([N+](=O)[O-])cc2)cc1'
    statistics_fragment_path= 'data/result/statistics_functional_1328/frag_smiles_main.xlsx'
    statistics_location_path='data/result/statistics_functional_1328/frag_location_main.xlsx'
    folder_path_a=highlight_chemical_atoms(folder_path,filename_without_extension, parent_dir,statistics_fragment_path,statistics_location_path)

