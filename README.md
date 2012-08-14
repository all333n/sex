## 简介
这个程序模拟N个二倍体个体，每个个体携带2个染色体并伴随着有害突变（有害突变的突变率为U，／代 ／染色体），并且有一定量的nbA位点受到性冲突选择。性别修饰基因位于染色体的正中。当一个个体繁殖时，减数分裂期间的染色体交换次数服从参数为L的泊松分布，并且交换事件平均的发生在整条染色体上。

模拟中首先要设定一个给定的预定代数使有性的比例达到稳定，然后使有性的比例开始进化。该程序会产生一个 _resultats.txt_ 文件，储存各参数的值和模拟花费的时间，另外一个文件—— _result\_N..txt_ 将会分别记录雌性和雄性的平均有性的比例，平均适合度每条染色体的平均有害突变数，性冲突位点的遗传多样性，固定的突变数和群体中的雄性数量。

这个程序同样产生了两个文件 _female.txt_ , _male,txt_ 给出了雌性和雄性不同时间点的分布的 _sigma_ 值。

文件 _sexdiff.h_ 定义了函数的原型和全局变量，其中 *chr* 结构代表染色体。主函数**in main.cpp**,调用了读取参数值 *in parametres.txt* 的函数,并且把参数值写入 *resultats.txt* 文件，然后利用函数 *recursion* 进行模拟。文件 *parametres.txt* 中的参数必须以\* 开始，并且每行代表一个参数集合，模拟程序依次进行模拟。文件 *SelRec.cpp* 包含计算个体适合度和重组染色体的子函数。

## 来源
This simulate program written by Denis Roze in the paper named *Differential Selection Between The sexes And Selection For Sex* .
