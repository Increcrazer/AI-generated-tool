% 生成顺序输出的伪随机数，并保存在xx_raw.txt文件里，每次生成的伪随机数具有随机性，需要及时保存
%% SDV=211
N = 256;
%% 周期S(偏振随机)
num_each = 64;

% 创建字符数组
random_string = [repmat('1', 1, num_each), repmat('3', 1, num_each), repmat('5', 1, num_each), repmat('7', 1, num_each)];

% 随机打乱字符顺序
random_string = random_string(randperm(N));

% 新建并打开.txt文件
fileID = fopen('00_raw.txt', 'w');
fprintf(fileID, '%s\n', random_string);
fclose(fileID);

%% 周期S+ S- SR SL
str = repmat('1', 1, 256);
fileID = fopen('01_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

str = repmat('3', 1, 256);
fileID = fopen('0C_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

str = repmat('5', 1, 256);
fileID = fopen('0D_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

str = repmat('7', 1, 256);
fileID = fopen('0E_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% 周期V+
str = repmat('C', 1, 256);
fileID = fopen('14_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% SV随机
characters = ['1','C'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('03_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% SDV随机
characters = ['1','1','8','C'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('04_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% SP+SN随机 SR+SL随机 SP+SN+SR+SL随机
characters = ['1','3'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('05_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

characters = ['5','7'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('0F_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

characters = ['1','3','5','7'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('06_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% 全随机
characters = ['1','3','5','7','8','9','A','B','C','D','E','F'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('00_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

fileID = fopen('02_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% S随机匹配P/N/R/L
characters = ['1','D','E','F'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('07_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

characters = ['3','C','E','F'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('08_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

characters = ['5','C','D','F'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('09_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

characters = ['7','C','D','E'];
% 计算集合中字符的数量
numChars = length(characters);
% 生成256个随机索引 (每个索引的概率相同)
randomIndices = randi([1, numChars], 1, 256);
% 根据索引从字符集合中选择对应的字符
str = characters(randomIndices);
fileID = fopen('0A_raw.txt', 'w');
fprintf(fileID, '%s\n', str);
fclose(fileID);

%% 30DP/DN/DR/DL+226VP
% 计算各个字符的数量
num_DP = 30;
% num_DN = 30;
% num_DR = 30;
% num_DL = 30;
num_VP = 226;

% 创建字符数组
random_string = repmat('C', 1, N);

% 30DP+256VP
random_string(1:num_DP) = '8';
random_string(num_DP+1:num_DP+num_VP) = 'C';
fileID = fopen('10_raw.txt', 'w');
fprintf(fileID, '%s', random_string);
fclose(fileID);

% 30DN+256VP
random_string(1:num_DP) = '9';
random_string(num_DP+1:num_DP+num_VP) = 'C';
fileID = fopen('11_raw.txt', 'w');
fprintf(fileID, '%s', random_string);
fclose(fileID);

% 30DR+256VP
random_string(1:num_DP) = 'A';
random_string(num_DP+1:num_DP+num_VP) = 'C';
fileID = fopen('12_raw.txt', 'w');
fprintf(fileID, '%s', random_string);
fclose(fileID);

% 30DL+256VP
random_string(1:num_DP) = 'B';
random_string(num_DP+1:num_DP+num_VP) = 'C';
fileID = fopen('13_raw.txt', 'w');
fprintf(fileID, '%s', random_string);
fclose(fileID);


