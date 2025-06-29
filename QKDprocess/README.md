
## statelabel_patternanalysis_SDV说明
- 首先会对脉冲强度进行统计，当arrset维度是3时说明统计正确，不正确的话需要修改count_resol
![](https://s2.loli.net/2025/06/09/piB8WAr6EMvIu5F.jpg)

- 然后会根据arrange_list进行搜索（该算法复制了PNRL搜索，实际上并不需要搜索，但笔者懒得改直接复用了）
![](https://s2.loli.net/2025/06/09/iWTyJzcIfbBke89.jpg)

- 代码还会输出模式效应的统计分析表格，如下所示：
![](https://s2.loli.net/2025/06/09/peytYqcOg8dRLIk.png)
关于模式效应的分析见:
Yuanfei Gao and Zhiliang Yuan, "Suppression of patterning effect using IQ modulator for high-speed quantum key distribution systems," Opt. Lett. 48, 1068-1071 (2023)
