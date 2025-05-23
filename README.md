# 复杂代码输入输出格式说明
## SNSPD_period_angle.m
表格格式：

![image](https://github.com/user-attachments/assets/208e22f9-c292-49eb-9bd1-f99b29b0a7ef)
## SIP_xxx.m
这些文件建议一起使用，用于生成prbs_generator需要的伪随机数码和每个GTX通道对应的伪随机数码
### SIP_randomgene
生成的伪随机数码格式如下：

![image](https://github.com/user-attachments/assets/303a702a-d8c5-4793-8767-98f9984a42e1)
### SIP_listprocess
将SIP_randomgene生成的伪随机数码转化为能直接输入给prbs_generator的格式。转化的伪随机数码格式如下：

![image](https://github.com/user-attachments/assets/a5409d3d-b79b-4046-9f27-e117fc7beb53)
### SIP_randomgene
根据“编码关系.xlsx”中的码字映射关系，将xx_raw.txt中的源码转换成各个GTX通道的输出，“编码关系.xlsx”格式如下：

![image](https://github.com/user-attachments/assets/38bfde6d-6d3f-43f4-bed5-fb7b7fd28fe2)

输出格式如下：

![image](https://github.com/user-attachments/assets/c791ec6b-51db-48de-943f-8977f116ea69)


