from rdkit import Chem
from rdkit.Chem import Draw

def annotate_atoms(smiles):
    # 解析 SMILES
    mol = Chem.MolFromSmiles(smiles)

    # 标注原子编号
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx()+1)
    # 生成带有原子编号的 SMILES
    annotated_smiles = Chem.MolToSmiles(mol)
        # 生成化学式图片
    img = Draw.MolToImage(mol, legend=annotated_smiles)

    return annotated_smiles, img



if __name__ == "__main__":


    #smiles="COC(C=C1)=CC2=C1C(C=CN=C3C(O)=O)=C3C2.CC4=C(C=O)C=CC=C4"
    smiles="COC(C=C1)=CC2=C1C(C=CN=C3C(O)=O)=C3C2.CC4=C(C([H])=O)C=CC=C4"

    # 获取带有标注的 SMILES 和图片
    annotated_smiles, img = annotate_atoms(smiles)
    # 打印标注的 SMILES
    print("Annotated SMILES:", annotated_smiles)
    # 显示图片
    img.show()