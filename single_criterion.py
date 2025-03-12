from remove_arrary_outerside import remove_outer_layer
import matplotlib.pyplot as plt
from binarize_by_ratio import binarize_by_ratio
from binary_functional_highatom import remove_rows_and_columns
from syn_1_img import syn_1_jpeg
import numpy as np
import os
from tokenizer import tokenize_smiles
from highlight_num_atoms import highlight_num_atoms
from remove_no_alpha import find_and_remove_non_letters
from exact_smile_fragment import extract_specified_atoms
from find_functional_group import identify_functional_groups
from rdkit import Chem
import pandas as pd
from save_funcional_group import save_functional_excel
from openpyxl import load_workbook

def attn_grey_binary(attn,folder_name_img,mol_name,ratio,i,j):
    # input one mol 8*8 attn,
    # output img of the mol attn,jpg is attn ,png is binarize

    attn_single=attn[i][j]
    attn_del2=remove_outer_layer(attn_single)
    plt.imshow(attn_del2, cmap='Greys')
    mol_name_ticks=tokenize_smiles(mol_name)
    plt.xticks(range(len(mol_name_ticks)),mol_name_ticks)
    plt.yticks(range(len(mol_name_ticks)),mol_name_ticks)
    plt.xlabel(f'{i}_{j}')
    folder_name_mol=os.path.join(folder_name_img,mol_name)
    if not os.path.exists(folder_name_mol):
        os.makedirs(folder_name_mol)
    single_name_greys = f'{i}_{j}.jpg'
    file_path_greys =os.path.join(folder_name_mol, single_name_greys)
    plt.savefig(file_path_greys)
    plt.close()

    ## binarize attn array
    attn_del2_bin = binarize_by_ratio(attn_del2,ratio)  #ratio=0.95 # set front ratio to  0
    #attn_del2_bin = find_and_keep_largest_block(attn_del2_bin) #find max bolk

    plt.imshow(attn_del2_bin, cmap='Greys')
    plt.xticks(range(len(mol_name_ticks)),mol_name_ticks)
    plt.yticks(range(len(mol_name_ticks)),mol_name_ticks)
    plt.xlabel(f'{i}_{j}')
    single_name_binary = f'{i}_{j}.png'
    file_path_binary =os.path.join(folder_name_mol,single_name_binary)
    plt.savefig(file_path_binary)
    plt.close()

    single_name_npy=f'{i}_{j}.npy'
    file_path_single_npy = os.path.join(folder_name_mol,single_name_npy)
    np.save(file_path_single_npy, attn_del2_bin)
    return folder_name_mol,file_path_single_npy,single_name_npy

def highlight_chemical_atoms (folder_path,filename_without_extension,file_path_single_npy,single_name_npy):
    if 'H' in filename_without_extension:
        pass
    else:
        file_path=file_path_single_npy
        attn_binary = np.load(file_path)
        np.fill_diagonal(attn_binary,0)#set diagonal element to 0
        attn_binary_del=attn_binary

        ## del non-character
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
        # find no_zero_element
        non_zero_positions = np.nonzero(attn_binary_letter)
        non_zero_col=np.unique(non_zero_positions[0]).tolist()
        non_zero_row=np.unique(non_zero_positions[1]).tolist()
        # print("non_zero_col不为零元素的行索引：", non_zero_col)
        # print("non_zero_row不为零元素的列索引：", non_zero_row)

        ## find functonal group
        smiles=filename_without_extension
        atom_indices=non_zero_col
        frag_smiles= extract_specified_atoms(smiles, atom_indices)
        detected_groups = identify_functional_groups(frag_smiles,smiles)

        #set first  row and title for next all row
        # df = pd.DataFrame(detected_groups)
        # df.to_excel("functional_groups_main.xlsx", index=False)
        # print("官能团检测结果已保存到 functional_groups_main.xlsx")

        ##set all row
        df_excel=pd.read_excel('functional_groups_single.xlsx')
        df_new=pd.DataFrame(detected_groups)
        df_excel.loc[len(df_excel)] = df_new.iloc[0,:]
        df_excel.to_excel('functional_groups_single.xlsx', index=False)

        ## highlight atoms
        atom_indices = non_zero_row
        img_data = highlight_num_atoms(smiles, atom_indices)

        # 保存图像为文件
        suffix_h='highatom'
        file_name_no_suffix=os.path.splitext(single_name_npy)[0]
        img_filename=f'{file_name_no_suffix}.jpeg'
        img_path=os.path.join(folder_path,img_filename)
        with open(img_path, "wb") as f:
            f.write(img_data)
        print(file_name_no_suffix)
        # #显示图像
        # from PIL import Image
        # import io
        # img = Image.open(io.BytesIO(img_data))
        # img.show()
        return folder_path,atom_indices

if __name__ == "__main__":
    ratio=0.95
    i=0
    j=0
    folder_name_criterion='result/attention/first/criterion'
    attn = np.load('./result/attention/first/npy/NC(=S)Nc1ccccc1.npy')
    filename_without_extension = 'NC(=S)Nc1ccccc1'
    filename_without_extension = filename_without_extension.replace("x", "/")
    filename_without_extension_re = tokenize_smiles(filename_without_extension)  # re syn like cl to one element
    filename_without_extension_forsaving = filename_without_extension.replace("/", "x")  # this filename for saving address
    ## show

    img_64_filename,file_path_single_npy, single_name_npy = attn_grey_binary(attn, folder_name_criterion, filename_without_extension,ratio,i,j)
    syn_name_attn = f'{filename_without_extension_forsaving}_{ratio}.jpg'
    output_image_attn = os.path.join(folder_name_criterion, syn_name_attn)
    attn_syn_img = syn_1_jpg(img_64_filename, syn_name_attn)  # syn 64 attention to 1 image
    print(img_64_filename)

    suffix_b = '_binarize'
    syn_name_binarize = f'{filename_without_extension_forsaving}{suffix_b}_{ratio}.jpg'
    output_image_binarize = os.path.join(folder_name_criterion, syn_name_binarize)
    attn_syn_binarize = syn_1_png(img_64_filename, syn_name_binarize)

    annotated_smiles, img = annotate_atoms(filename_without_extension)
    suffix_a = '_mol'
    syn_name_mol = f'{filename_without_extension_forsaving}{suffix_a}_{ratio}.jpg'
    output_image_mol = os.path.join(img_64_filename,syn_name_mol)
    img.save(output_image_mol)

    suffix_c = '_lightatom'
    folder_path, atom_indices= highlight_chemical_atoms(img_64_filename, filename_without_extension,file_path_single_npy,single_name_npy)
    syn_name_atom = f'{filename_without_extension_forsaving}{suffix_c}_{ratio}.jpg'
    output_image_atom = os.path.join(folder_name_criterion, syn_name_atom)
    attn_syn_atom = syn_1_jpeg(img_64_filename, syn_name_atom)

