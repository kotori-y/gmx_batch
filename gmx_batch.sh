#!/bin/bash

# 使用set -e确保脚本遇到错误时退出
set -e

# 设置错误捕捉，发生错误时打印错误信息并继续执行
trap 'echo "Error occurred at line $LINENO"' ERR

# 使用sobtop生成小分子（*.mol2）所需的gro，top，以及itp文件
# sobtop下载地址：http://sobereva.com/soft/Sobtop/
preprocess_mol() {
  mol_file_path=$1
  mol_file_dir=$2
  mol_name=$3

  tgt_dir="$mol_file_dir/../kotori_results/$mol_name/"
  if [ ! -d "$tgt_dir" ]; then
    mkdir -p "$tgt_dir"
    echo "$tgt_dir has been created"
  fi

  gro_path="$tgt_dir/molecule.gro"
  top_path="$tgt_dir/molecule.top"
  itp_path="$tgt_dir/molecule.itp"

  # 使用sobtop
  printf "$mol_file_path\n2\n$gro_path\n1\n2\n4\n$top_path\n$itp_path\n0\n" | ./sobtop

  # 在.ipt文件末尾加入以下3行
  echo '#ifdef POSERS' >> "$itp_path"
  echo '#include "posre.itp"' >> "$itp_path"
  echo '#endif' >> "$itp_path"
}

# 合并小分子和蛋白的gro文件，并将蛋白下的topol.top修改至小分子目录下
combine() {
  work_path=$1
  mol_name=$2

  mol_gro_path="$work_path/kotori_results/$mol_name/molecule.gro"
  protein_gro_path="$work_path/protein/protein.gro"

  n_mol_line=$(wc -l < "$mol_gro_path")  # 获取小分子文件行数
  n_mol_atom=$(sed -n "2p" "$mol_gro_path")  # 获取小分子原子总数
  n_protein_atom=$(sed -n "2p" "$protein_gro_path")  # 获取蛋白原子总数

  total_atom=$((n_protein_atom + n_mol_atom))  # 总原子数
  out_file="$work_path/kotori_results/$mol_name/complex.gro"

  # 写入文件内容
  head -n 1 "$protein_gro_path" > "$out_file"  # 写入蛋白文件的第1行
  echo "$total_atom" >> "$out_file"  # 写入总原子数
  head -n -1 "$protein_gro_path" | tail -n +3 >> "$out_file"  # 写入蛋白文件的第1至第倒数第3行
  head -n -1 "$mol_gro_path" | tail -n +3 >> "$out_file"  # 写入小分子文件的第1至第倒数第3行
  tail -n 1 "$protein_gro_path" >> "$out_file"  # 写入蛋白文件的最后一行

  src_topol_file="$work_path/protein/topol.top"
  tgt_topol_file="$work_path/kotori_results/$mol_name/topol.top"

  # 将蛋白引用改成相对小分子的路径
  sed '/; Include chain topologies/{n; s|^#include "|#include "../../protein/|; :a; n; /^; Include /!{ s|^#include "|#include "../../protein/|; ba; }}' "$src_topol_file" > "$tgt_topol_file"

  # 获取topol_Protein_chain_A.itp的行号，在这行之前添加小分子的引用
  line_number=$(grep -n '#include "topol_Protein_chain_A.itp"' "$src_topol_file" | cut -d ":" -f1)
  sed -i "${line_number}i #include \"molecule.itp\"" "$tgt_topol_file"

  # 添加小分子信息
  echo "$mol_name     1" >> "$tgt_topol_file"
}

# 生成position restraint文件
gen_posre_file() {
  work_path=$1
  mol_name=$2

  mol_gro_path="$work_path/kotori_results/$mol_name/molecule.gro"
  out_file="$work_path/kotori_results/$mol_name/posre.itp"
  printf "0\n" | gmx genrestr -f "$mol_gro_path" -o "$out_file"
}

# 运行模拟，捕捉任何错误
run() {
  work_path=$1
  mol_name=$2
  gpu=$3

  mol_root="$work_path/kotori_results/$mol_name"

  box_in="$mol_root/complex.gro"
  box_out="$mol_root/protein_box.gro"

  water_out="$mol_root/protein_sol.gro"
  topol_file="$mol_root/topol.top"

  tpr_out="$mol_root/em.tpr"
  iron_out="$mol_root/system.gro"

  em_gro="$mol_root/em.gro"
  eq_tpr="$mol_root/eq.tpr"
  eq_gro="$mol_root/eq.gro"
  md_tpr="$mol_root/md.tpr"

  # 捕捉错误的关键点
  {
    # 加盒子
    gmx editconf -f "$box_in" -o "$box_out" -d 1.0 -bt cubic && \
    # 溶剂化
    gmx solvate -cp "$box_out" -o "$water_out" -p "$topol_file" && \
    # 生成初步结构文件，需要忽略一个warning，使用时请注意
    gmx grompp -f "$work_path/em.mdp" -c "$water_out" -r "$water_out" -p "$topol_file" -o "$tpr_out" -maxwarn 1 && \
    # 离子化，选择15，即SOL
    printf "15\n" | gmx genion -s "$tpr_out" -p "$topol_file" -o "$iron_out" -neutral -conc 0.15 && \
    # 能量最小化准备
    gmx grompp -f "$work_path/em.mdp" -c "$iron_out" -r "$iron_out" -p "$topol_file" -o "$tpr_out" && \
    # 能量最小化
    pushd "$mol_root" && gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 16 -gpu_id "$gpu" && popd && \
    # 预平衡准备
    gmx grompp -f "$work_path/eq.mdp" -c "$em_gro" -p "$topol_file" -o "$eq_tpr" -r "$em_gro" -maxwarn 1 && \
    # 预平衡
    pushd "$mol_root" && gmx mdrun -v -deffnm eq -ntmpi 1 -ntomp 16 -gpu_id "$gpu" && popd && \
    # 动力学
    gmx grompp -f "$work_path/md.mdp" -c "$eq_gro" -r "$eq_gro" -p "$topol_file" -o "$md_tpr" && \
    pushd "$mol_root" && gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 16 -gpu_id "$gpu" && popd && echo 'done'
  } || {
    echo "An error occurred during the run step for $mol_name. Skipping to the next molecule."
    return 1  # 确保报错时返回，并不影响下一次循环
  }
}

# 获得分析结果
analyse{
  work_path=$1
  mol_name=$2

  mol_root="$work_path/kotori_results/$mol_name"

  md_trr="$mol_root/md.trr"
  md_xtc="$mol_root/md.xtc"
  md_edr="$mol_root/md.edr"
  md_tpr="$mol_root/md.tpr"

# 提取xtc轨迹(xtc精度稍低，缺少每个原子的速度和受力信息)
gmx trjconv -f $md_trr -o $md_xtc && \
# 提取体系total能量
printf "15\n" | gmx energy -f $md_edr -o "$mol_root/${mol_name}_energy.xvg" && \

# 计算Protein rmsd,选两次protein
printf "1\n1\n" | gmx rms -f $md_xtc -s $md_tpr -o "$mol_root/${mol_name}_rmsd_protein.xvg" && \
# 计算骨架 rmsd,选两次backbone
printf "4\n4\n" | gmx rms -f $md_xtc -s $md_tpr -o "$mol_root/${mol_name}_rmsd_backbone.xvg" && \
# 计算Mol-pro rmsd，第一次protein，第二次mol，代表结构对蛋白做叠合，但只计算mol的rmsd
printf "1\n13\n" | gmx rms -f $md_xtc -s $md_tpr -o "$mol_root/${mol_name}_rmsd_lig.xvg" && \
# 计算Mol rmsd，第一次mol，第二次mol
printf "13\n13\n" | gmx rms -f $md_xtc -s $md_tpr -o "$mol_root/${mol_name}_rmsd_lig_lig.xvg" && \
# 氢键数目，第一次pro，第二次mol
printf "1\n13\n" | gmx hbond -f $md_xtc -s $md_tpr -num "$mol_root/${mol_name}_hbnum.xvg" && \
}

# 主函数
main() {
  work_path=$1
  mol_file_path=$2
  gpu_id=$3

  mol_file_name=$(basename "$mol_file_path")
  mol_file_dir=$(dirname "$mol_file_path")
  mol_name="${mol_file_name%.*}"

  preprocess_mol "$mol_file_path" "$mol_file_dir" "$mol_name" && \
    combine "$work_path" "$mol_name" && \
    gen_posre_file "$work_path" "$mol_name" && \
    run "$work_path" "$mol_name" "$gpu_id" && \
    analyse "$work_path" "$mol_name"
}

# 批量处理，确保每次循环执行完毕后重新初始化
work_path=$1
mol2_dir=$2
gpu_id=$3

find "$work_path/$mol2_dir" -type f | while read -r file; do
    main "$work_path" "$file" "$gpu_id"
done
