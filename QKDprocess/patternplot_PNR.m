%% 初始化设置 - PR Applied风格
clear; clc;
close all;

% 文件设置
file_pattern = 'SDV-*-ch12_result_order1.xlsx';  % 修改为PNR前缀
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
    
    % 确定目标列（使用新的PNR分组函数）
    target_columns = generate_columns_for_order_pnr(current_order);
    
    % 保留存在的列
    target_columns = target_columns(ismember(target_columns, data.Properties.VariableNames));
    data = data(:, target_columns);
    
    % 归一化处理 - 按最后一位分组分别归一化
    data = normalize_by_last_char(data, target_columns);
    
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
    target_columns = generate_columns_for_order_pnr(current_order);
    target_columns = target_columns(ismember(target_columns, order_data.Properties.VariableNames));
    order_data = order_data(:, target_columns);
    
    if height(order_data) == 0
        continue;
    end
    
    %% 按最后一位字母分组处理
    last_chars = cellfun(@(x) x(end), target_columns, 'UniformOutput', false);
    unique_last_chars = unique(last_chars);
    
    %% 创建每个组的单独图形
    for group_idx = 1:length(unique_last_chars)
        last_char = unique_last_chars{group_idx};
        
        % 获取该组的所有列
        group_columns = target_columns(strcmp(last_chars, last_char));
        group_data = order_data(:, group_columns);
        
        % 跳过没有数据的组
        if isempty(group_columns)
            continue;
        end
        
        % 为该组创建颜色方案
        num_patterns = length(group_columns);
        group_colors = lines(num_patterns);  % 使用lines颜色方案
        
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
        
        for i = 1:length(group_columns)
            col_name = group_columns{i};
            col_data = group_data{:, col_name};
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
                        'Color', group_colors(i, :));

                line_handles(end+1) = h;
                legend_entries{end+1} = sprintf('%s (n=%d)', col_name, length(valid_data));
            end
        end

        %% 图表美化 - PR风格
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
        title(ax, sprintf('Pattern-Dependent Photon Statistics (Order %d, %c-ending)', current_order, last_char), ...
              'FontSize', pr_fontsize+2, ...
              'FontWeight', 'bold', ...
              'FontName', pr_fontname);
        
        % 设置合理的X轴范围
        all_valid_data = group_data{:, :};
        all_valid_data = all_valid_data(all_valid_data > 0 & ~isnan(all_valid_data));
        if ~isempty(all_valid_data)
            xlim(ax, [min(all_valid_data)*0.9, max(all_valid_data)*1.1]);
        end
        
        % 添加基准线
        plot(ax, [1 1], ylim(ax), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        
        % 网格
        grid(ax, 'on');
        grid(ax, 'minor');
        ax.MinorGridLineStyle = ':';
        ax.MinorGridAlpha = 0.2;
        
        % 图例 - 放在左上角
        if ~isempty(line_handles)
            legend(ax, line_handles, legend_entries, ...
                   'Location', 'northwest', ...
                   'FontSize', pr_fontsize-1, ...
                   'FontName', pr_fontname, ...
                   'Box', 'off');
        end
        
        %% 保存主图
        png_filename = fullfile(output_dir, sprintf('Fig_Order%d_%c_ending.png', current_order, last_char));
        saveas(fig, png_filename);
        
        % 保存高分辨率版本
        pdf_filename = fullfile(output_dir, sprintf('Fig_Order%d_%c_ending.pdf', current_order, last_char));
        exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
        
        fprintf('  Saved: %s\n', png_filename);
        fprintf('  Saved: %s (vector)\n', pdf_filename);
        
        %% 保存数据表格和统计信息
        writetable(group_data, fullfile(data_output_dir, sprintf('Order%d_%c_ending_Data.xlsx', current_order, last_char)));
        writetable(stats_table, fullfile(data_output_dir, sprintf('Order%d_%c_ending_Statistics.xlsx', current_order, last_char)));
        
        fprintf('  Group %c-ending: Processed %d patterns with total %d pulses\n', ...
                last_char, length(group_columns), sum(stats_table.Count));
        
    end
end

%% 保存合并数据和汇总信息
all_data_file = fullfile(data_output_dir, 'All_Orders_Data.xlsx');
writetable(all_norm_data, all_data_file);
fprintf('\nSaved combined data: %s\n', all_data_file);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Figures saved to: %s\n', output_dir);
fprintf('Data saved to: %s\n', data_output_dir);

%% 辅助函数
function columns = generate_columns_for_order_pnr(order)
    state_types = {'P', 'N', 'R'};
    
    if order == 0
        % 0阶：所有单态
        columns = {'P', 'N', 'R'};
        columns = sort_columns_by_pnr(columns);
        return;
    end
    
    % 生成所有模式
    all_patterns = generate_combinations_pnr(state_types, order + 1);
    
    % 按最后一位分组并排序
    columns = group_and_sort_patterns_by_last_char(all_patterns);
end

function combinations = generate_combinations_pnr(elements, n)
    if n == 0
        combinations = {''};
    elseif n == 1
        combinations = elements;
    else
        sub_combinations = generate_combinations_pnr(elements, n-1);
        combinations = {};
        for i = 1:length(elements)
            for j = 1:length(sub_combinations)
                combinations{end+1} = [elements{i} sub_combinations{j}];
            end
        end
    end
    combinations = combinations(:);
end

function sorted_columns = sort_columns_by_pnr(columns)
    pattern_chars = char(columns);
    num_patterns = zeros(length(columns), size(pattern_chars, 2));
    
    for i = 1:length(columns)
        for j = 1:size(pattern_chars, 2)
            switch pattern_chars(i, j)
                case 'P'
                    num_patterns(i, j) = 1;
                case 'N'
                    num_patterns(i, j) = 2;
                case 'R'
                    num_patterns(i, j) = 3;
            end
        end
    end
    
    [~, idx] = sortrows(num_patterns);
    sorted_columns = columns(idx);
end

function sorted_patterns = group_and_sort_patterns_by_last_char(patterns)
    % 按最后一位分组
    last_chars = cellfun(@(x) x(end), patterns, 'UniformOutput', false);
    unique_last_chars = unique(last_chars);
    
    sorted_patterns = {};
    
    % 对每个分组单独排序
    for i = 1:length(unique_last_chars)
        last_char = unique_last_chars{i};
        
        % 找出该组的所有模式
        group_patterns = patterns(strcmp(last_chars, last_char));
        
        % 按PNR顺序排序该组模式
        group_patterns = sort_columns_by_pnr(group_patterns);
        
        % 添加到结果中
        sorted_patterns = [sorted_patterns; group_patterns(:)];
    end
end

function normalized_data = normalize_by_last_char(data, columns)
    % 按最后一位分组进行归一化
    last_chars = cellfun(@(x) x(end), columns, 'UniformOutput', false);
    unique_last_chars = unique(last_chars);
    
    normalized_data = data;
    
    for i = 1:length(unique_last_chars)
        last_char = unique_last_chars{i};
        
        % 获取该组的所有列
        group_columns = columns(strcmp(last_chars, last_char));
        
        % 找出该组的全同模式作为基准
        baseline_pattern = repmat(last_char, 1, length(group_columns{1}));
        baseline_col = find(strcmp(group_columns, baseline_pattern));
        
        if ~isempty(baseline_col)
            baseline_data = data{:, group_columns{baseline_col}};
            valid_baseline = baseline_data(baseline_data > 0 & ~isnan(baseline_data));
            
            if ~isempty(valid_baseline)
                current_baseline_mean = mean(valid_baseline);
                
                % 归一化该组的所有列
                for col_idx = 1:length(group_columns)
                    col_name = group_columns{col_idx};
                    col_data = data{:, col_name};
                    normalized_data{:, col_name} = col_data / current_baseline_mean;
                end
            end
        end
    end
end