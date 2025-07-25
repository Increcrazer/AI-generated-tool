## delay_accuracy_analysis说明
首先得有示波器存储的波形文件，代码中是'1250M_PM.csv'

![alt text](waveform.jpg)

然后会放大某个码字对应的波形，并用红色虚线标出波形平坦区：

![alt text](partialwaveform.jpg)

然后会放大平坦区，并对平坦区做平滑，添加三种不同位置的高斯光脉冲：

![alt text](eopartialwaveform.jpg)

然后给出三种位置高斯光脉冲经过调制后的误码率，以此来说明延时精度的必要性：

![alt text](result.jpg)

## statelabel_patternanalysis_SDV说明
首先会对脉冲强度进行统计，当arrset维度是3时说明统计正确，不正确的话需要修改count_resol

![](https://s2.loli.net/2025/06/09/piB8WAr6EMvIu5F.jpg)

然后会根据arrange_list进行搜索（该算法复制了PNRL搜索，实际上并不需要搜索，但笔者懒得改直接复用了）

![](https://s2.loli.net/2025/06/09/iWTyJzcIfbBke89.jpg)

代码还会输出模式效应的统计分析表格，如下所示：

![](https://s2.loli.net/2025/06/09/peytYqcOg8dRLIk.png)

关于模式效应的分析见:
Yuanfei Gao and Zhiliang Yuan, "Suppression of patterning effect using IQ modulator for high-speed quantum key distribution systems," Opt. Lett. 48, 1068-1071 (2023)
