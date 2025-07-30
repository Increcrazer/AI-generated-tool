## delay_accuracy_analysis说明
首先得有示波器存储的波形文件，代码中是'1250M_PM.csv'

![](https://s2.loli.net/2025/07/31/Y6Qrtls8xHeAZG3.png)

然后会放大某个码字对应的波形，并用红色虚线标出波形平坦区：

![](https://s2.loli.net/2025/07/31/zSQp4nbYHuGwD9M.png)

然后会放大平坦区，并对平坦区做平滑，添加三种不同位置的高斯光脉冲：

![](https://s2.loli.net/2025/07/31/zrqO15GUKywHcgS.png)

然后给出三种位置高斯光脉冲经过调制后的误码率，以此来说明延时精度的必要性：

![](https://s2.loli.net/2025/07/31/9J3tnWwfQMquDVG.png)

## statelabel_patternanalysis_SDV说明
首先会对脉冲强度进行统计，当arrset维度是3时说明统计正确，不正确的话需要修改count_resol

![](https://s2.loli.net/2025/06/09/piB8WAr6EMvIu5F.jpg)

然后会根据arrange_list进行搜索（该算法复制了PNRL搜索，实际上并不需要搜索，但笔者懒得改直接复用了）

![](https://s2.loli.net/2025/06/09/iWTyJzcIfbBke89.jpg)

代码还会输出模式效应的统计分析表格，如下所示：

![](https://s2.loli.net/2025/06/09/peytYqcOg8dRLIk.png)

关于模式效应的分析见:
Yuanfei Gao and Zhiliang Yuan, "Suppression of patterning effect using IQ modulator for high-speed quantum key distribution systems," Opt. Lett. 48, 1068-1071 (2023)

## statelabel_patternanalysis_SDV_patch说明
将statelabel_patternanalysis_SDV封装为一个函数，但是不输出任何图像，可以批量输出模式效应的统计分析表格，效果如下：

![](https://s2.loli.net/2025/07/31/YemhBsPGaQE2g31.png)

## patternplot_SDV说明
对statelabel_patternanalysis_SDV_patch输出的一些列表格进行进一步处理，以将这些大批量表格先各自按照SS的平均值计数进行归一化，再进行拼接，这样可以避免长时间采数导致的功率漂移也被算入pattern effect的影响，然后再绘制pattern频数图：

![](https://s2.loli.net/2025/07/31/dNDyjvs87A5IeWG.png)



