把以下三个对接复合物 PDB 放到本目录并保持命名一致:

  input/alllhrc_complex.pdb
  input/fllhttr_complex.pdb
  input/ylsllqr_complex.pdb

它们是不同肽与 4ey6 (乙酰胆碱酯酶) 在 PAS 位点的最优对接构象。
脚本 run_all.sh 会按前缀(<前缀>_complex.pdb)读取并逐体系运行 MD。
