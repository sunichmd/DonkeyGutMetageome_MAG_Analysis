# 配对
# -s 是否标准化 -t 过滤低丰度标准 -N 计算z、c score -M 计算模块 -S 点的大小范围
script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata.txt -n Group \
-p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0.005 -M FALSE -N F -w 14 -e 7 -o .0.005

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata1.txt -n Group \
-l taxonomy.txt -L 6 -p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0.005 -M FALSE -N T -w 14 -e 7 -o .0.005.taxonomy

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata1.txt -n Group \
-p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0 -M FALSE -N F -w 14 -e 7 -S 2,2 -o .all

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata1.txt -n Group \
-l taxonomy.txt -L 6 -p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0 -M FALSE -N T -w 14 -e 7 -S 2,2 -o .all.taxonomy

# 全部
# -s 是否标准化 -t 过滤低丰度标准 -N 计算z、c score -M 计算模块 -S 点的大小范围
script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata.txt -n Group \
-p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0.005 -M FALSE -N F -w 14 -e 7 -o .0.005

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata.txt -n Group \
-l taxonomy.txt -L 6 -p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0.005 -M FALSE -N T -w 14 -e 7 -o .0.005.taxonomy

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata.txt -n Group \
-p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0 -M FALSE -N F -w 14 -e 7 -S 2,2 -o .all

script/network_analyse_ansi.R -i MAGs_relative_abundance_abundance.xls -d metadata.txt -n Group \
-l taxonomy.txt -L 6 -p 0.05 -r 0.85 -m spearman -a Bonferroni -s F \
-t 0 -M FALSE -N T -w 14 -e 7 -S 2,2 -o .all.taxonomy