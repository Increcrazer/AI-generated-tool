%% 初始化设置 - PR Applied风格
clear; clc;
close all;

% 文件设置
file_pattern = 'SDV-8192random30s*_result_order0.xlsx';
output_dir = 'Processed_Figures';
data_output_dir = 'Processed_Data';

% 创建输出文件夹
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end
if ~exist(data_output_dir, 'dir')
    mkdir(data_output_dir); 
end

% 设置风格参数
pr_linewidth = 1.5;
pr_markersize = 6;
pr_fontsize = 11;
pr_fontname = 'Helvetica'; % PR期刊推荐字体

%% 智能文件处理函数
process_excel = @(filename) ...
    readtable(filename, ...
             'Sheet', 1, ...
             'VariableNamingRule', 'preserve', ...
             'NumHeaderLines', 0);

%% 批量处理文件
file_list = dir(file_pattern);
all_norm_data = table();

fprintf('=== Processing files ===\n');

for i = 1:length(file_list)
    filename = file_list(i).name;
    fprintf('Processing: %s\n', filename);
    
    % 提取order信息
    token = regexp(filename, 'order(\d+)', 'tokens');
    if ~isempty(token)
        current_order = str2double(token{1}{1});
    else
        current_order = 1;
    end
    
    % 读取数据
    data = process_excel(filename);
    
    % 移除统计行
    valid_rows = height(data) - 5;
    data = data(1:valid_rows, :);
    
    % 确定目标列
    if current_order == 1
        target_columns = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};
    else
        target_columns = generate_columns_for_order(current_order);
    end
    
    % 保留存在的列
    target_columns = target_columns(ismember(target_columns, data.Properties.VariableNames));
    data = data(:, target_columns);
    
    % 归一化处理
    baseline_column = repmat('S', 1, current_order + 1);
    if ismember(baseline_column, target_columns)
        baseline_data = data.(baseline_column);
        valid_baseline = baseline_data(baseline_data > 0 & ~isnan(baseline_data));
        
        if ~isempty(valid_baseline)
            current_baseline_mean = mean(valid_baseline);
            
            % 归一化所有列
            for col_idx = 1:length(target_columns)
                col_name = target_columns{col_idx};
                col_data = data.(col_name);
                data.(col_name) = col_data / current_baseline_mean;
            end
        end
    end
    
    % 添加元数据
    data.File = repmat({filename}, height(data), 1);
    data.Order = repmat(current_order, height(data), 1);
    
    % 合并数据
    all_norm_data = [all_norm_data; data];
end

%% 按order分组处理和可视化
unique_orders = unique(all_norm_data.Order);

for order_idx = 1:length(unique_orders)
    current_order = unique_orders(order_idx);
    
    % 筛选当前order数据
    order_mask = all_norm_data.Order == current_order;
    order_data = all_norm_data(order_mask, :);
    
    % 移除元数据列
    order_data(:, {'File', 'Order'}) = [];
    
    % 确定列顺序
    if current_order == 1
        target_columns = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};
    else
        target_columns = generate_columns_for_order(current_order);
    end
    
    target_columns = target_columns(ismember(target_columns, order_data.Properties.VariableNames));
    order_data = order_data(:, target_columns);
    
    if height(order_data) == 0
        continue;
    end
    
    %% 颜色方案（风格）
    num_patterns = length(target_columns);
    pr_colors = turbo(num_patterns);  % 自动生成足够多的颜色

    
    %% 创建风格的主图
    fig = figure('Position', [100, 100, 800, 500], ...
                 'Color', 'white', ...
                 'PaperPositionMode', 'auto');
    
    % 设置坐标轴
    ax = axes('Parent', fig, ...
              'Box', 'on', ...
              'LineWidth', 1.0, ...
              'FontSize', pr_fontsize, ...
              'FontName', pr_fontname, ...
              'TickDir', 'out', ...
              'TickLength', [0.015, 0.015]);
    hold(ax, 'on');
    
    %% 绘制分布图
    bin_width = 0.001; % 适当增加bin宽度以获得更平滑的分布
    legend_entries = {};
    line_handles = [];
    
    % 统计信息存储
    stats_table = table();
    
    for i = 1:length(target_columns)
        col_name = target_columns{i};
        col_data = order_data{:, col_name};
        valid_data = col_data(col_data > 0 & ~isnan(col_data));

        if ~isempty(valid_data)
            % 计算统计信息
            stats_table.ColName{i} = col_name;
            stats_table.Mean(i) = mean(valid_data);
            stats_table.Std(i) = std(valid_data);
            stats_table.Count(i) = length(valid_data);
            stats_table.FWHM(i) = 2 * sqrt(2 * log(2)) * std(valid_data); % 近似FWHM

            % 计算直方图（使用计数模式，不进行归一化）
            [counts, edges] = histcounts(valid_data, 'BinWidth', bin_width);
            centers = edges(1:end-1) + diff(edges)/2;

            % 应用高斯滤波平滑
            smoothed_counts = smoothdata(counts, 'gaussian', 5);

            % 绘制平滑后的计数曲线
            h = plot(ax, centers, smoothed_counts, ...
                    'LineWidth', pr_linewidth, ...
                    'Color', pr_colors(mod(i-1, size(pr_colors,1))+1, :));

            line_handles(end+1) = h;
            legend_entries{end+1} = sprintf('%s (n=%d)', col_name, length(valid_data));
        end
    end

    %% 图表美化 - 风格
    % 坐标轴标签
    xlabel(ax, 'Normalized Photon Counts', ...
           'FontSize', pr_fontsize+1, ...
           'FontWeight', 'bold', ...
           'FontName', pr_fontname);
    ylabel(ax, 'Pulse Counts', ...
           'FontSize', pr_fontsize+1, ...
           'FontWeight', 'bold', ...
           'FontName', pr_fontname);

    % 标题
    title(ax, sprintf('Pattern-Dependent Photon Statistics (Order %d)', current_order), ...
          'FontSize', pr_fontsize+2, ...
          'FontWeight', 'bold', ...
          'FontName', pr_fontname);
    
    % 网格
    grid(ax, 'on');
    grid(ax, 'minor');
    ax.MinorGridLineStyle = ':';
    ax.MinorGridAlpha = 0.2;
    
    % 图例 - 放在左上角
    if ~isempty(line_handles)
        legend(ax, line_handles, legend_entries, ...
               'Location', 'northwest', ...  % 修改为左上角
               'FontSize', pr_fontsize-1, ...
               'FontName', pr_fontname, ...
               'Box', 'off');
    end
    
    %% 保存主图（只保存PNG格式）
    % 高质量PNG
    png_filename = fullfile(output_dir, sprintf('Fig_Order%d_Distribution.png', current_order));
    saveas(fig, png_filename);
    fprintf('  Saved: %s\n', png_filename);
    
    %% 保存数据表格
    writetable(order_data, fullfile(data_output_dir, sprintf('Order%d_Data.xlsx', current_order)));
    
    % 保存统计信息
    writetable(stats_table, fullfile(data_output_dir, sprintf('Order%d_Statistics.xlsx', current_order)));
    
    fprintf('Order %d: Processed %d patterns with total %d pulses\n', ...
            current_order, length(target_columns), sum(stats_table.Count));
end

%% 保存合并数据和汇总信息
all_data_file = fullfile(data_output_dir, 'All_Orders_Data.xlsx');
writetable(all_norm_data, all_data_file);
fprintf('\nSaved combined data: %s\n', all_data_file);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Figures saved to: %s\n', output_dir);
fprintf('Data saved to: %s\n', data_output_dir);

%% 辅助函数
function columns = generate_columns_for_order(order)
    state_types = {'S', 'D', 'V'};
    columns = {};
    
    all_prefixes = generate_combinations(state_types, order);
    
    for p = 1:length(all_prefixes)
        prefix = all_prefixes{p};
        columns{end+1} = [prefix 'S'];
        columns{end+1} = [prefix 'D'];
    end
    
    columns = sort_columns_by_sdv(columns);
end

function combinations = generate_combinations(elements, n)
    % 处理n=0的情况
    if n == 0
        combinations = {''};  % 空字符串表示没有前缀
        return;
    elseif n == 1
        combinations = elements;
    else
        sub_combinations = generate_combinations(elements, n-1);
        combinations = {};
        for i = 1:length(elements)
            for j = 1:length(sub_combinations)
                combinations{end+1} = [elements{i} sub_combinations{j}];
            end
        end
    end
    combinations = combinations(:);
end

function sorted_columns = sort_columns_by_sdv(columns)
    pattern_chars = char(columns);
    num_patterns = zeros(length(columns), size(pattern_chars, 2));
    
    for i = 1:length(columns)
        for j = 1:size(pattern_chars, 2)
            switch pattern_chars(i, j)
                case 'S'
                    num_patterns(i, j) = 1;
                case 'D'
                    num_patterns(i, j) = 2;
                case 'V'
                    num_patterns(i, j) = 3;
            end
        end
    end
    
    [~, idx] = sortrows(num_patterns);
    sorted_columns = columns(idx);
end
