% 清除工作区和关闭所有图形窗口
clear all;
close all;

% 获取当前文件夹路径
currentFolder = pwd;

% 查找所有 MPD[纯数字].mat 文件
fileList = dir(fullfile(currentFolder, 'MPD*.mat'));
numFiles = length(fileList);

% 初始化有效文件列表
validFiles = {};

% 检查文件名是否严格匹配 MPD[纯数字].mat
for i = 1:numFiles
    fileName = fileList(i).name;
    % 使用正则表达式匹配 MPD + 纯数字 + .mat
    if ~isempty(regexp(fileName, '^MPD\d+\.mat$', 'once'))
        validFiles{end+1} = fileName;
    end
end

numValidFiles = length(validFiles);
k = zeros(1, numValidFiles);
D = zeros(1, numValidFiles);

% 检查是否有有效文件
if numValidFiles == 0
    error('未找到符合命名规则的 MPD[数字].mat 文件。');
end

% 创建新图形窗口
figure;
hold on;
grid on;
xlabel('Power (mW)');
ylabel('Current (uA)');
title('MPD');

% 循环处理每个有效文件
for i = 1:numValidFiles
    fileName = validFiles{i};
    filePath = fullfile(currentFolder, fileName);
    
    % 提取纯数字部分（例如 MPD22.mat → 22）
    mpdNum = str2double(regexp(fileName, '\d+', 'match', 'once'));
    varName = ['MPD' num2str(mpdNum)]; % 变量名，如 MPD1、MPD22
    
    % 加载数据
    data = load(filePath);
    vars = fieldnames(data);
    mpdData = data.(vars{1});

    % 将数据存入工作区（可选）
    assignin('base', varName, mpdData);
    
    % 提取第一行（x 轴）和第二行（y 轴）
    xData = mpdData(1, :);
    yData = mpdData(2, :);
    
    % 计算斜率 k 和线性度偏差 D
    k(i) = (yData(end) - yData(1)) / (xData(end) - xData(1));
    D(i) = calculate_linearity_deviation(mpdData);
    
    % 绘制曲线，并在图例中显示 k 和 D
    legendEntry = sprintf('%s (k=%.2f, D=%.4f)', varName, k(i), D(i));
    plot(xData, yData, 'DisplayName', legendEntry, 'LineWidth', 1.5);
end

% 添加图例并调整位置
legend('show', 'Location', 'best');


function D = calculate_linearity_deviation(data)
    % 提取第一行和第二行数据
    row1 = 10*log10(data(1, :));    % mW转换为dBm，为了和GD文档统一
    row2 = data(2, :);
    
    % 找到第一行中最接近0.8和0.2的值及其列索引
    [~, c1] = min(abs(row1 - 10*log10(0.8)));
    [~, c2] = min(abs(row1 - 10*log10(0.2)));
    
    % 获取对应的P1, P2, I1, I2
    P1 = row1(c1);
    P2 = row1(c2);
    I1 = row2(c1);
    I2 = row2(c2);
    
    % 计算delta和D
    delta = (P1 - P2) - 10 * log10(I1 / I2);
    D = 10^(delta / 10) - 1;
end

% function D = calculate_linearity_deviation(data)  % 等效写法
%     % 提取第一行和第二行数据
%     row1 = data(1, :);    
%     row2 = data(2, :);
%     
%     % 找到第一行中最接近0.8和0.2的值及其列索引
%     [~, c1] = min(abs(row1 - 0.8));
%     [~, c2] = min(abs(row1 - 0.2));
%     
%     % 获取对应的P1, P2, I1, I2
%     P1 = row1(c1);
%     P2 = row1(c2);
%     I1 = row2(c1);
%     I2 = row2(c2);
%     
%     D = P1/P2/(I1/I2)-1;
% end

