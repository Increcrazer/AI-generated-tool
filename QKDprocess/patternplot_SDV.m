%% 初始化设置
clear; clc;
file_pattern = 'SDV-8192random30s*_result_order2.xlsx'; % 根据实际情况调整
output_dir = 'Smart_Processed_Results2_random1'; % 输出文件夹

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

% 用于存储每个文件的order信息
file_orders = zeros(length(file_list), 1);

for i = 1:length(file_list)
    filename = file_list(i).name;
    fprintf('正在智能处理: %s\n', filename);
    
    % 1. 从文件名提取order信息
    % 假设文件名格式: *_result_orderX.xlsx
    token = regexp(filename, 'order(\d+)', 'tokens');
    if ~isempty(token)
        current_order = str2double(token{1}{1});
    else
        % 如果文件名中没有order信息，使用默认值
        current_order = 1;
        warning('文件名 %s 中未找到order信息，使用默认order=1', filename);
    end
    file_orders(i) = current_order;
    
    % 2. 读取原始数据（跳过第1行标题）
    data = process_excel(filename);
    
    % 3. 自动识别统计行（最后5行）
    valid_rows = height(data) - 5; % 总行数-5
    data = data(1:valid_rows, :);  % 保留有效数据
    
    % 4. 根据order确定目标列
    % 对于order=1: SS, DS, VS, SD, DD, VD
    % 对于order=2: SSS, DSS, VSS, SDS, DDS, VDS, SSV, DSV, VSV, SDV, DDV, VDV, SSD, DSD, VSD, SDD, DDD, VDD
    % 等等...
    
    if current_order == 1
        target_columns = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};
    else
        % 生成对应order的所有模式列名
        target_columns = generate_columns_for_order(current_order);
    end
    
    % 5. 检查数据中是否存在这些列
    missing_columns = setdiff(target_columns, data.Properties.VariableNames);
    if ~isempty(missing_columns)
        warning('文件 %s 缺少以下列: %s', filename, strjoin(missing_columns, ', '));
        target_columns = intersect(target_columns, data.Properties.VariableNames);
    end
    
    % 6. 提取目标列
    data = data(:, target_columns);
    
    % 7. 对当前文件的数据进行归一化
    % 确定基准列（全S模式）
    baseline_column = repmat('S', 1, current_order + 1);
    if ismember(baseline_column, target_columns)
        current_baseline_mean = mean(data.(baseline_column)(data.(baseline_column) > 0 & ~isnan(data.(baseline_column))));
    else
        % 如果没有基准列，使用第一个S结尾的列
        s_ending_columns = target_columns(cellfun(@(x) x(end) == 'S', target_columns));
        if ~isempty(s_ending_columns)
            current_baseline_mean = mean(data.(s_ending_columns{1})(data.(s_ending_columns{1}) > 0 & ~isnan(data.(s_ending_columns{1}))));
        else
            current_baseline_mean = 1;
            warning('文件 %s 未找到合适的基准列进行归一化', filename);
        end
    end
    
    % 归一化当前文件的数据
    norm_data = data;
    for col = target_columns
        col_name = col{1};
        norm_data.(col_name) = norm_data.(col_name) / current_baseline_mean;
    end
    
    % 8. 添加文件标识列
    norm_data.File = repmat({filename}, height(norm_data), 1);
    norm_data.Order = repmat(current_order, height(norm_data), 1);
    
    % 9. 合并归一化后的数据
    all_norm_data = [all_norm_data; norm_data];
end

%% 按order分组处理
unique_orders = unique(file_orders);

for order_idx = 1:length(unique_orders)
    current_order = unique_orders(order_idx);
    
    % 筛选当前order的数据
    order_mask = all_norm_data.Order == current_order;
    order_data = all_norm_data(order_mask, :);
    
    % 移除标识列
    order_data(:, {'File', 'Order'}) = [];
    
    % 确定当前order的目标列
    target_columns = generate_columns_for_order(current_order);
    order_data = order_data(:, intersect(target_columns, order_data.Properties.VariableNames));
    
    % 检查是否有数据
    if height(order_data) == 0
        warning('Order=%d 没有有效数据，跳过', current_order);
        continue;
    end
    
    % 生成图例标签
    legend_labels = order_data.Properties.VariableNames;
    
    % 设置颜色
    colors = lines(width(order_data));
    
    %% 脉冲计数直方图
    % 参数设置
    bin_width = 0.0001; % 保持与之前相同的bin宽度
    normalization = 'count'; % 改为计数模式
    
    figure('Position', [100, 100, 800, 600], 'Name', sprintf('Order %d', current_order));
    hold on;
    
    % 预分配图形对象
    bar_handles = gobjects(1, width(order_data)); 
    
    % 计算总脉冲数（用于后续百分比显示）
    total_pulses = sum(sum(~isnan(order_data{:,1:width(order_data)})));
    
    for i = 1:width(order_data)
        col_name = order_data.Properties.VariableNames{i};
        col_data = order_data{:,i};
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
                                   'DisplayName', sprintf('%s (n=%d)', col_name, length(valid_data)));
        end
    end
    
    % 图表美化
    title(sprintf('Pulse Count Distribution by Pattern (Order=%d)', current_order), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Energy (a.u.)', 'FontSize', 12);
    ylabel('Pulse Count', 'FontSize', 12);
    
    % 优化图例显示
    legend('show', 'Location', 'northwest', 'FontSize', 10);
    grid on;
    box on;
    
    % 可选：添加总脉冲数标注
    text(0.8, 0.98, sprintf('Total Pulses: %d\nOrder: %d', total_pulses, current_order), ...
         'Units', 'normalized', ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 11, ...
         'BackgroundColor', [1 1 1 0.7]);
    
    %% 保存结果
    saveas(gcf, fullfile(output_dir, sprintf('Order%d_Normalized_Distribution.png', current_order)));
    
    % 保存当前order的数据
    writetable(order_data, fullfile(output_dir, sprintf('Order%d_Combined_Data.xlsx', current_order)));
end

% 保存所有数据的汇总
writetable(all_norm_data, fullfile(output_dir, 'All_Orders_Combined_Data.xlsx'));

fprintf('\n智能处理完成！结果已保存至: %s\n', output_dir);
fprintf('处理的Order值: %s\n', mat2str(unique_orders'));

%% 辅助函数：根据order生成列名
function columns = generate_columns_for_order(order)
    state_types = {'S', 'D', 'V'};
    columns = {};
    
    % 生成所有前order位的组合
    all_prefixes = generate_combinations(state_types, order);
    
    % 生成所有模式（最后一位只能是S或D）
    for p = 1:length(all_prefixes)
        prefix = all_prefixes{p};
        % 添加以S结尾的模式
        columns{end+1} = [prefix 'S'];
        % 添加以D结尾的模式
        columns{end+1} = [prefix 'D'];
    end
    
    % 按照SDV顺序排序
    columns = sort_columns_by_sdv(columns);
end

%% 辅助函数：生成所有组合
function combinations = generate_combinations(elements, n)
    if n == 1
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

%% 辅助函数：按SDV顺序排序列名
function sorted_columns = sort_columns_by_sdv(columns)
    % 将SDV映射为数字：S=1, D=2, V=3
    pattern_chars = char(columns);
    
    % 创建一个数值矩阵用于排序
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
    
    % 使用sortrows按数值矩阵排序
    [~, idx] = sortrows(num_patterns);
    sorted_columns = columns(idx);
end
