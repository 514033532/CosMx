# 加载所需库
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# 读取数据文件
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_type/proportion/Group1_DN2-VS-DN1/CellType_prop_DE_list_in_Group1_DN2-VS-DN1.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

读取文件
文件部分如下
	W	Pvalue	Significance	DN1_Mean	DN2_Mean	Mean_difference	Change	B1	B3	B4	B5	B6
Monocyte.cell	1	0.4	Not	5.455	6.55333333333333	1.09833333333333	up in case	4.18	7.36	3.84	8.12	7.07


分组信息
样本	分组
B1	DN2
B3	DN2
B4	DN1
B5	DN2
B6	DN1
F1	NC
F2	NC
F3	NC
F4	NC
F5	NC
F6	NC
F7	NC
A2	NC
D4	NC
D5	NC
D6	NC


比较"DN1", "DN2")两组的差异，分别独立绘制图形，标注P值，每个细胞独立保存PDF和png，去网格，保留XY，文章发表格式，XY，字体14，标题XY黑体20，只给出P小于0.05的，字体不考虑，且独立图是个封闭的框架图



