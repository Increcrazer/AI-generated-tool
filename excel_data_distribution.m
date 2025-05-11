function excel_data_distribution(filename, sheet)
    % ANALYZE_EXCEL_DATA 分析Excel文件中第一列数据的统计特性
    %
    % 功能描述：
    %   该函数读取指定Excel文件中指定工作表的第一列数据，
    %   计算其基本统计特性（包括最大值、最小值、均值、方差等），
    %   并绘制带有统计参数标注的分布直方图。
    %
    % 输入参数：
    %   filename - 字符串，Excel文件名（包含路径）
    %   sheet    - 字符串或数值，可选参数，指定工作表名称或索引（默认为1）
    %
    % 输出：
    %   在命令行窗口显示统计结果，并生成两个图形窗口：
    %   1. 带有统计参数标注的直方图
    %   2. 数据箱线图
    %
    % 示例：
    %   analyze_excel_data('data.xlsx')          % 分析data.xlsx的第一个工作表
    %   analyze_excel_data('data.xlsx', 'Sheet1') % 分析指定工作表
    %
    % 创建日期：2023-10-20
    % 作者：DeepSeek Chat
    
    % 参数检查与默认值设置
    if nargin < 2
        sheet = 1; % 默认使用第一个工作表
    end
    
    % 数据读取与预处理
    try
        % 读取Excel第一列数据
        data = readtable(filename, 'Sheet', sheet, 'Range', 'A:A');
        
        % 转换为数组并移除NaN值
        vector = table2array(data);
        vector = vector(~isnan(vector));
        
        if isempty(vector)
            error('第一列未找到有效的数值数据');
        end
        
        % 计算统计量
        data_min = min(vector);
        data_max = max(vector);
        data_mean = mean(vector);
        data_median = median(vector);
        data_variance = var(vector);
        data_std = std(vector);
        
        % 命令行输出统计结果
        fprintf('\n====== 数据统计结果 ======\n');
        fprintf('数据点数: %d\n', length(vector));
        fprintf('最小值: %.4f\n', data_min);
        fprintf('最大值: %.4f\n', data_max);
        fprintf('平均值: %.4f\n', data_mean);
        fprintf('中位数: %.4f\n', data_median);
        fprintf('方差: %.4f\n', data_variance);
        fprintf('标准差: %.4f\n', data_std);
        fprintf('==========================\n\n');
        
        % 创建带有统计参数的标题文本
        stats_text = sprintf(['统计参数:  N=%d\n'...
                            'Min=%.2f  Max=%.2f\n'...
                            'Mean=%.2f  Med=%.2f\n'...
                            'Var=%.2f  Std=%.2f'], ...
                            length(vector), data_min, data_max, ...
                            data_mean, data_median, ...
                            data_variance, data_std);
        
        % 绘制直方图（带统计参数标注）
        figure;
        histfit(vector); % 带正态分布拟合的直方图
        title({'数据分布直方图'; stats_text});
        xlabel('插损');
        ylabel('频数');
        grid on;
        
        % 绘制箱线图
        figure;
        boxplot(vector);
        title('数据箱线图');
        ylabel('插损');
        set(gca, 'FontSize', 10); % 设置字体大小
        
    catch ME
        % 错误处理
        fprintf('\n错误发生: %s\n', ME.message);
        fprintf('请检查:\n');
        fprintf('1. 文件路径是否正确\n');
        fprintf('2. 文件是否被其他程序占用\n');
        fprintf('3. 工作表名称或索引是否正确\n');
    end
end
