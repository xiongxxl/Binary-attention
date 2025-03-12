# from rdkit import Chem
# from rdkit.Chem import Draw
# # 定义分子
# smiles = 'CCO'  # 示例 SMILES 结构
# mol = Chem.MolFromSmiles(smiles)
# # 设置需要高亮的原子索引
# highlight_atoms = [0, 1]  # 例如高亮第0和第1个原子
# # 设置高亮颜色（蓝色：RGB形式）
# highlight_color = {atom_idx: (1.0, 0.0, 1.0) for atom_idx in highlight_atoms}
# # 绘制分子图，指定高亮原子及其颜色
# img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_color)
from rdkit import Chem
from rdkit.Chem import Draw
# 定义分子
smiles = 'CCO'  # 示例 SMILES 结构
mol = Chem.MolFromSmiles(smiles)
# 设置需要高亮的原子索引
highlight_atoms = [0, 1]  # 例如高亮第0和第1个原子
# 设置高亮颜色（蓝色：RGB形式）
highlight_color = {atom_idx: (0.0, 0.0, 1.0) for atom_idx in highlight_atoms}
# 创建 MolDraw2DCairo 对象
drawer = Draw.MolDraw2DCairo(300, 300)  # 300x300 像素大小的图像
# 开始绘制
drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_color)
# 完成并获取图像
drawer.FinishDrawing()
img = drawer.GetDrawingText()
# 将绘制的图像保存到文件
with open("molecule_highlighted.png", "wb") as f:
    f.write(img)