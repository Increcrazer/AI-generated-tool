% 设定文件前缀
prefix = 'SDV-8192random30s';
suffix = '.csv';

% 其他参数
bin_width = 16; % [ps]
freq = 1.25 *10^9;  % [Hz]
count_resol = 25;  % count resolution
order = 2; % 阶数（0: 单态, 1: 一阶, 2: 二阶, ...）
arrange_list = char("DSDDSSSVSSSSSSVVSDSSDVSVSD");

% 获取当前目录下所有符合的文件
files = dir([prefix, '*', suffix]);

% 依次处理每个文件
for k = 1:length(files)
    filename = files(k).name;
    disp(['处理文件: ', filename]);
    try
        analyze_pulse_patterns(filename, bin_width, freq, count_resol, arrange_list, order);
    catch ME
        warning(['处理失败: ', filename, ', 错误信息: ', ME.message]);
    end
end

disp('批量处理完成。');


%% 
function analyze_pulse_patterns(filename, bin_width, freq, count_resol, arrange_list, order)
    % 分析脉冲模式并生成统计结果
    % 输入参数:
    %   filename: CSV文件名 (字符串)
    %   bin_width: 时间宽度 [ps] (数值)
    %   freq: 频率 [Hz] (数值)
    %   count_resol: 计数分辨率 (数值)
    %   arrange_list: 状态序列 (字符数组或字符串)
    %   order: 阶数 (0: 单态, 1: 一阶, 2: 二阶, ...) (数值)

    %% find states of the pulses 
    period = 1 / freq * 10^12;    % [ps]
    MINPEAKDISTANCE = period / bin_width - 3;

    time = readmatrix(filename, 'Range', '1:1');
    data = readmatrix(filename, 'Range', '2:2');
    [~, index_list] = findpeaks(data, 'MINPEAKHEIGHT', 1, 'MINPEAKDISTANCE', MINPEAKDISTANCE);  
    index_list = index_list(2:end - 1);
    index_first = index_list(1);
    index_last = index_list(end);

    gate_ratio = 320 / 800;
    pulse = zeros(1, length(index_list));
    for i = 1:length(index_list)
        pulse(i) = sum(data(index_list(i) - period / bin_width / 2 * gate_ratio:1:index_list(i) + period / bin_width / 2 * gate_ratio));
    end

    log_pulse = log10(pulse);
    [y,x] = hist(log_pulse,count_resol);
    nonzero_index = find(y); % find 函数用于查找数组中非零元素的索引
    arrset = find_continuous_sequences(nonzero_index);
    if (numel(arrset) ~= 3)
        disp('请调整分辨率');
        return;
    end

    resolution = max(log_pulse) / count_resol;
    state1_range = [x(arrset{1}(1)) - resolution, x(arrset{1}(end)) + resolution];
    state2_range = [x(arrset{2}(1)) - resolution, x(arrset{2}(end)) + resolution];
    state3_range = [x(arrset{3}(1)) - resolution, x(arrset{3}(end)) + resolution];

    % code the pulses
    code_pulse = zeros(1, length(index_list));
    for i = 1:length(index_list)
        if (log_pulse(i) >= state1_range(1) && log_pulse(i) <= state1_range(2))
            code_pulse(i) = 1;
        elseif (log_pulse(i) >= state2_range(1) && log_pulse(i) <= state2_range(2))
            code_pulse(i) = 2;
        elseif (log_pulse(i) >= state3_range(1) && log_pulse(i) <= state3_range(2))
            code_pulse(i) = 3;
        else
            code_pulse(i) = 0;
        end
    end

    % find pulse position
    dic1 = struct('S', 1, 'D', 2, 'V', 3);
    dic2 = struct('S', 1, 'D', 3, 'V', 2);
    dic3 = struct('S', 2, 'D', 1, 'V', 3);
    dic4 = struct('S', 2, 'D', 3, 'V', 1);
    dic5 = struct('S', 3, 'D', 1, 'V', 2);
    dic6 = struct('S', 3, 'D', 2, 'V', 1);

    m = [];
    for i = 1:6
        for j = 1:length(arrange_list)
            m = [m, eval(strcat('dic', num2str(i))).(arrange_list(j))];
        end
        eval([strcat('a', num2str(i)), '= m;']);
        m = [];
    end

    arrange_list_order = 0;
    selected_dic = 0;
    for i = 1:6
        m = eval(strcat('a', num2str(i)));
        for j = 1:length(code_pulse) - length(arrange_list)
            pulse_slice = code_pulse(j:j + length(arrange_list) - 1);
            pulse_slice_index = find(code_pulse(j:j + length(arrange_list) - 1));
            if (m(pulse_slice_index) == pulse_slice(pulse_slice_index))
                disp('Find position!');
                disp(strcat('Dictionary is dic', num2str(i)));
                disp(strcat('a(1) corresponds to code_pulse(', num2str(j), ')'));
                arrange_list_order = j;
                selected_dic = i;
                break;
            else
                continue;
            end
        end
        if arrange_list_order > 0
            break;
        end
    end

    %% 生成所有可能的模式组合并按正确顺序排序
    state_types = {'S', 'D', 'V'};
    pattern_count = 3^order * 2; % 前order位可以是任意状态，最后一位只能是S或D

    % 首先生成所有前order位的组合
    all_prefixes = generate_combinations(state_types, order);

    % 生成所有模式
    all_patterns = {};
    for p = 1:length(all_prefixes)
        prefix = all_prefixes{p};
        % 添加以S结尾的模式
        all_patterns{end + 1} = [prefix 'S'];
        % 添加以D结尾的模式
        all_patterns{end + 1} = [prefix 'D'];
    end

    % 按照特殊规则排序：
    % 1. 首先按最后一位分组（S结尾的在前，D结尾的在后）
    % 2. 在每个分组内，按前缀的字典顺序排序
    pattern_names = cell(1, pattern_count);
    s_patterns = {}; % S结尾的模式
    d_patterns = {}; % D结尾的模式

    % 分离S结尾和D结尾的模式
    for p = 1:length(all_patterns)
        pattern = all_patterns{p};
        if pattern(end) == 'S'
            s_patterns{end + 1} = pattern;
        else
            d_patterns{end + 1} = pattern;
        end
    end

    % 对S结尾的模式按SDV顺序排序
    s_patterns = sort_patterns_by_sdv(s_patterns);

    % 对D结尾的模式按SDV顺序排序  
    d_patterns = sort_patterns_by_sdv(d_patterns);

    % 合并：先放所有S结尾的，再放所有D结尾的
    idx = 1;
    for p = 1:length(s_patterns)
        pattern_names{idx} = s_patterns{p};
        idx = idx + 1;
    end

    for p = 1:length(d_patterns)
        pattern_names{idx} = d_patterns{p};
        idx = idx + 1;
    end

    fprintf('阶数: %d, 模式数量: %d\n', order, pattern_count);

    %% 将pulse序列归类并写入Excel，计算完整统计量
    % 创建输出矩阵（增加统计行）
    output_data = cell(length(code_pulse) + 5, 3 + pattern_count);

    % 列标题
    headers = {'Pulse Number', 'Pulse Counts', 'Pulse Type'};
    for j = 1:pattern_count
        headers{end + 1} = pattern_names{j};
    end

    % 使用找到的dic{i}进行归类
    current_dic = eval(['dic' num2str(selected_dic)]);

    % 填充前三列数据
    for pulse_idx = 1:length(code_pulse)
        output_data{pulse_idx, 1} = pulse_idx;
        output_data{pulse_idx, 2} = pulse(pulse_idx);

        current_code = code_pulse(pulse_idx);
        if current_code == current_dic.S
            output_data{pulse_idx, 3} = 'S';
        elseif current_code == current_dic.D
            output_data{pulse_idx, 3} = 'D';
        elseif current_code == current_dic.V
            output_data{pulse_idx, 3} = 'V';
        else
            output_data{pulse_idx, 3} = 'Unknown';
        end
    end

    % 初始化数据收集
    pattern_counts = cell(pattern_count, 1);

    % 检测模式并收集数据
    for pulse_idx = 1:length(code_pulse) - order
        % 获取当前窗口的状态序列
        window_states = cell(1, order + 1);
        window_codes = zeros(1, order + 1);

        for k = 0:order
            state_code = code_pulse(pulse_idx + k);
            window_codes(k + 1) = state_code;
            if state_code == current_dic.S
                window_states{k + 1} = 'S';
            elseif state_code == current_dic.D
                window_states{k + 1} = 'D';
            elseif state_code == current_dic.V
                window_states{k + 1} = 'V';
            end
        end

        % 检查窗口是否包含有效状态
        if all(window_codes > 0)
            % 生成模式名称
            pattern_name = '';
            for k = 1:order + 1
                pattern_name = [pattern_name window_states{k}];
            end

            % 查找模式在列表中的位置
            pattern_idx = find(strcmp(pattern_names, pattern_name));

            if ~isempty(pattern_idx)
                last_count = pulse(pulse_idx + order);
                pattern_counts{pattern_idx} = [pattern_counts{pattern_idx}; last_count];
                output_data{pulse_idx + order, pattern_idx + 3} = last_count;
            end
        end
    end

    % 计算统计量
    stats_row = length(code_pulse) + 1;

    % 初始化统计数组
    means = zeros(1, pattern_count);
    stds = zeros(1, pattern_count);
    valid_patterns = ~cellfun(@isempty, pattern_counts);

    % 计算基本统计量
    for p = 1:pattern_count
        if valid_patterns(p)
            means(p) = mean(pattern_counts{p});
            stds(p) = std(pattern_counts{p});
        end
    end

    % 计算每个模式的归一化因子和偏差
    norm_means = zeros(1, pattern_count);
    norm_stds = zeros(1, pattern_count);
    deviations = zeros(1, pattern_count);

    % 查找基准模式
    all_s_pattern = '';  % 全S模式，用于S结尾的模式
    all_d_pattern = '';  % 全D模式，用于D结尾的模式
    for k = 1:order + 1
        all_s_pattern = [all_s_pattern 'S'];
        all_d_pattern = [all_d_pattern 'D'];
    end

    s_baseline_idx = find(strcmp(pattern_names, all_s_pattern));
    d_baseline_idx = find(strcmp(pattern_names, all_d_pattern));

    % 修正：使用if-else语句
    if ~isempty(s_baseline_idx) && valid_patterns(s_baseline_idx)
        s_baseline_mean = means(s_baseline_idx);
    else
        s_baseline_mean = 1;
    end

    if ~isempty(d_baseline_idx) && valid_patterns(d_baseline_idx)
        d_baseline_mean = means(d_baseline_idx);
    else
        d_baseline_mean = 1;
    end

    for p = 1:pattern_count
        if valid_patterns(p)
            pattern_name = pattern_names{p};

            % 根据最后一位选择基准
            if pattern_name(end) == 'S'
                baseline_mean = s_baseline_mean;
                baseline_name = all_s_pattern;
            else % pattern_name(end) == 'D'
                baseline_mean = d_baseline_mean;
                baseline_name = all_d_pattern;
            end

            % 计算归一化值和偏差
            norm_means(p) = means(p) / baseline_mean;
            norm_stds(p) = stds(p) / baseline_mean;
            deviations(p) = (means(p) - baseline_mean) / baseline_mean * 100;
        end
    end

    % 填充统计行
    output_data{stats_row, 1} = 'Mean (counts)';
    output_data{stats_row + 1, 1} = 'Std Dev (counts)';
    output_data{stats_row + 2, 1} = 'Normalized Mean';
    output_data{stats_row + 3, 1} = 'Normalized Std Dev';
    output_data{stats_row + 4, 1} = 'Deviation Percentage';

    for p = 1:pattern_count
        if valid_patterns(p)
            output_data{stats_row, p + 3} = means(p);
            output_data{stats_row + 1, p + 3} = stds(p);
            output_data{stats_row + 2, p + 3} = norm_means(p);
            output_data{stats_row + 3, p + 3} = norm_stds(p);
            output_data{stats_row + 4, p + 3} = sprintf('%.3f%%', deviations(p));
        end
    end

    % 写入Excel文件
    output_filename = strrep(filename, '.csv', sprintf('_result_order%d.xlsx', order));

    % 将标题和数据合并
    all_data = [headers; output_data];

    % 写入Excel
    try
        writecell(all_data, output_filename, 'Sheet', 1);
        fprintf('数据已成功保存到: %s\n', output_filename);
    catch ME
        error('写入Excel失败: %s\n请确保文件没有被其他程序打开。', ME.message);
    end

    fprintf('\n结果已写入文件: %s\n', output_filename);

end

%% 连续数字序列查找函数
function arrset = find_continuous_sequences(nonzero_index)
    arrset = cell(0, 0);
    if isempty(nonzero_index)
        return;
    end

    start_idx = 1;
    n = numel(nonzero_index);

    while start_idx <= n
        end_idx = start_idx;
        while (end_idx < n) && (nonzero_index(end_idx) + 1 == nonzero_index(end_idx + 1))
            end_idx = end_idx + 1;
        end

        arrset{end + 1} = nonzero_index(start_idx:end_idx);
        start_idx = end_idx + 1;
    end
end

%% 生成所有组合函数（递归）
function combinations = generate_combinations(elements, n)
    if n == 1
        combinations = elements;
    else
        sub_combinations = generate_combinations(elements, n - 1);
        combinations = {};
        for i = 1:length(elements)
            for j = 1:length(sub_combinations)
                combinations{end + 1} = [elements{i} sub_combinations{j}];
            end
        end
    end
    combinations = combinations(:);
end

%% 排序函数
function sorted_patterns = sort_patterns_by_sdv(patterns)
    % 将SDV映射为数字：S=1, D=2, V=3
    pattern_chars = char(patterns);
    
    % 创建一个数值矩阵用于排序
    num_patterns = zeros(length(patterns), size(pattern_chars, 2));
    
    for i = 1:length(patterns)
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
    sorted_patterns = patterns(idx);
end
