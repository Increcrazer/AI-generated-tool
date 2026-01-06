%% 初始化设置
clear; clc;
file_pattern = 'SDV-8192random30s*_result.xlsx'; % 根据实际情况调整
output_dir = 'Smart_Processed_Results_random1'; % 输出文件夹
target_columns = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'}; % 目标列（SS/DS/VS/SD/DD/VD）

% 创建输出文件夹
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end

%% 智能文件处理函数
process_excel = @(filename) ...
    readtable(filename, ...
             'Sheet', 1, ...
             'VariableNamingRule', 'preserve', ...
             'NumHeaderLines', 0); % 不用跳过标题

%% 批量处理文件
file_list = dir(file_pattern);
all_norm_data = table();

for i = 1:length(file_list)
    filename = file_list(i).name;
    fprintf('正在智能处理: %s\n', filename);
    
    % 1. 读取原始数据（跳过第1行标题）
    data = process_excel(filename);
    
    % 2. 自动识别统计行（最后5行）
    valid_rows = height(data) - 5; % 总行数-5
    data = data(1:valid_rows, :);  % 保留有效数据
    
    % 3. 提取目标列
    data = data(:, target_columns);
    
    % 4. 对当前文件的数据进行归一化
    % 计算当前文件的SS列非零均值
    current_ss_mean = mean(data.SS(data.SS > 0 & ~isnan(data.SS)));
    
    % 归一化当前文件的数据
    norm_data = data;
    for col = target_columns
        col_name = col{1};
        norm_data.(col_name) = norm_data.(col_name) / current_ss_mean;
    end
    
    % 5. 合并归一化后的数据
    all_norm_data = [all_norm_data; norm_data];
end

%% 第一张图：S态和D态的脉冲数量直方图
% 参数设置
bin_width = 0.0001;
normalization = 'count'; % 改为计数模式

figure('Position', [100, 100, 800, 600], 'Name', 'S-state and D-state Count Distribution');
hold on;

% 定义S态和D态列
s_state_columns = {'SS', 'DS', 'VS'};
d_state_columns = {'SD', 'DD', 'VD'};

% 颜色设置
colors = lines(2); % 只需要2种颜色：S态和D态
s_color = colors(1, :);
d_color = colors(2, :);

% 计算S态数据（合并SS、DS、VS）
s_state_data = [];
for i = 1:length(s_state_columns)
    col_data = all_norm_data{:, s_state_columns{i}};
    valid_data = col_data(col_data > 0 & ~isnan(col_data));
    s_state_data = [s_state_data; valid_data];
end

% 计算D态数据（合并SD、DD、VD）
d_state_data = [];
for i = 1:length(d_state_columns)
    col_data = all_norm_data{:, d_state_columns{i}};
    valid_data = col_data(col_data > 0 & ~isnan(col_data));
    d_state_data = [d_state_data; valid_data];
end

% 预分配图形对象
bar_handles_state = gobjects(1, 2);

% 绘制S态分布
if ~isempty(s_state_data)
    [counts_s, edges_s] = histcounts(s_state_data, 'BinWidth', bin_width, 'Normalization', normalization);
    centers_s = edges_s(1:end-1) + diff(edges_s)/2;
    bar_handles_state(1) = stairs(centers_s, counts_s, ...
                                 'Color', s_color, ...
                                 'LineWidth', 2.5, ...
                                 'DisplayName', sprintf('S-state (n=%d)', length(s_state_data)));
end

% 绘制D态分布
if ~isempty(d_state_data)
    [counts_d, edges_d] = histcounts(d_state_data, 'BinWidth', bin_width, 'Normalization', normalization);
    centers_d = edges_d(1:end-1) + diff(edges_d)/2;
    bar_handles_state(2) = stairs(centers_d, counts_d, ...
                                 'Color', d_color, ...
                                 'LineWidth', 2.5, ...
                                 'DisplayName', sprintf('D-state (n=%d)', length(d_state_data)));
end

% 图表美化
title('Pulse Count Distribution: S-state vs D-state', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Energy (a.u.)', 'FontSize', 12);
ylabel('Pulse Count', 'FontSize', 12);
legend('show', 'Location', 'northwest', 'FontSize', 10);
grid on;
box on;

% 计算并显示统计信息
s_state_count = length(s_state_data);
d_state_count = length(d_state_data);
total_state_count = s_state_count + d_state_count;
text(0.98, 0.98, sprintf('Total: %d\nS-state: %d (%.1f%%)\nD-state: %d (%.1f%%)', ...
     total_state_count, s_state_count, s_state_count/total_state_count*100, ...
     d_state_count, d_state_count/total_state_count*100), ...
     'Units', 'normalized', ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 10, ...
     'BackgroundColor', [1 1 1 0.7]);

% 保存第一张图
saveas(gcf, fullfile(output_dir, 'S_and_D_state_Count_Distribution.png'));

%% 第二张图：所有通道的脉冲计数分布
% 参数设置
bin_width = 0.0001; % 保持与之前相同的bin宽度
normalization = 'count'; % 计数模式

figure('Position', [100, 100, 800, 600], 'Name', 'All Channels Count Distribution');
hold on;

colors = lines(width(all_norm_data));
legend_labels = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};

% 预分配图形对象
bar_handles = gobjects(1, width(all_norm_data)); 

% 计算总脉冲数（用于后续百分比显示）
total_pulses = sum(sum(~isnan(all_norm_data{:,1:width(all_norm_data)})));

for i = 1:width(all_norm_data)
    col_data = all_norm_data{:,i};
    valid_data = col_data(col_data > 0 & ~isnan(col_data));
    
    if ~isempty(valid_data)
        % 计算直方图（使用计数模式）
        [counts, edges] = histcounts(valid_data, 'BinWidth', bin_width, 'Normalization', normalization);
        
        % 转换为阶梯图中心点
        centers = edges(1:end-1) + diff(edges)/2;
        
        % 绘制条形图（阶梯状）
        bar_handles(i) = stairs(centers, counts, ...
                               'Color', colors(i,:), ...
                               'LineWidth', 2.5, ...
                               'DisplayName', sprintf('%s (n=%d)', legend_labels{i}, length(valid_data)));
    end
end

% 图表美化
title('Pulse Count Distribution by first-order Pattern', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Energy (a.u.)', 'FontSize', 12);
ylabel('Pulse Count', 'FontSize', 12);

% 优化图例显示
legend('show', 'Location', 'northwest', 'FontSize', 10);
grid on;
box on;

% 可选：添加总脉冲数标注
text(0.8, 0.98, sprintf('Total Pulses: %d', total_pulses), ...
     'Units', 'normalized', ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 11, ...
     'BackgroundColor', [1 1 1 0.7]);

% 保存第二张图
saveas(gcf, fullfile(output_dir, 'Smart_Normalized_Distribution.png'));

%% 保存数据
writetable(all_norm_data, fullfile(output_dir, 'Smart_Combined_Data.xlsx'));

fprintf('智能处理完成！结果已保存至: %s\n', output_dir);
