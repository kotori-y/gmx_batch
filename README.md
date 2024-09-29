# GROMACS 批处理

## 环境

- 系统：Linux
- Sobtop：1.0

## 使用方法

### 下载Sobtop

从[这里](http://sobereva.com/soft/Sobtop/)获取Sobtop，然后将本脚本置于Sobtop目录。在Sobtop目录下运行：

```shell
chmod +x *
```

### 准备文件

1. 准备一个工作目录，首先放入含有蛋白相关文件的文件夹，并将其命名为<code>protein</code>，下面是一个示例：

   ```plaintext 
   work_path/
   └── protein/
       ├── posre_Protein_chain_A.itp
       ├── posre_Protein_chain_B.itp
       ├── posre_Protein_chain_C.itp
       ├── posre_Protein_chain_D.itp
       ├── posre_Protein_chain_E.itp
       ├── topol_Protein_chain_A.itp
       ├── topol_Protein_chain_B.itp
       ├── topol_Protein_chain_C.itp
       ├── topol_Protein_chain_D.itp
       ├── topol_Protein_chain_E.itp
       ├── protein.gro
       └── topol.top
   ```

2. 将<code>gmx_batch.sh</code>第57行的<code>#include "topol_Protein_chain_A.itp"</code>替换成<code>topol.top</code>文件中 <code>; Include chain topologies</code>下面一行的内容。

3. 接着，在工作目录下放入只含有待计算小分子的文件夹，命名不限（下面假设为<code>mol2</code>），分子格式为<code>.mol2</code>。

   ```plaintext 
   work_path/
   ├── protein/
   │   ...
   │   └── topol.top
   └── mol2/
       ├── Molecule1.mol2
       ...
       └── Molecule100.mol2
   ```

   **注意：分子的文件名应与<code>.mol2</code>文件一致**

3. 在工作目录下放入<code>em.mdp</code>，<code>eq.mdp</code>以及<code>md.mdp</code>:

   ```plaintext 
   work_path/
   ├── protein/
   │   ...
   │   └── topol.top
   ├── mol2/
   │   ...
   │   └── Molecule100.mol2
   ├── em.mdp
   ├── eq.mdp
   └── md.mdp
   ```

### 运行

```shell
./gmx_batch.sh /path/to/workpath $mol_dir_name $gpu_id
```

其中，<code>/path/to/workpath</code>为工作目录相对Sobtop目录的路径，<code>\$mol_dir_name</code>为小分子文件夹名，<code>\$gpu_id</code>为显卡ID。

<del>什么，你说你没有GPU？那请自行修改脚本中gpu的相关参数</del>

## TODO

- [ ] 蛋白文件批量生成，实现仅需<code>.pdb</code>与<code>.mol2</code>，一行命令即可运行完整分子动力学。

- [ ] 开箱即用的可交互GROMACS软件。


## Kotori reminds you

  ことりは世界で一番かわいい