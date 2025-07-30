% 设定文件前缀
prefix = 'SDV-1024random1-40s';
suffix = '.csv';

% 定义伪随机序列
arrange_list = char("SDSVVDVSDDDSSSSVVSDSVVSVSDVDVVD");

% 获取当前目录下所有符合的文件
files = dir([prefix, '*', suffix]);

% 依次处理每个文件
for k = 1:length(files)
    filename = files(k).name;
    disp(['处理文件: ', filename]);
    try
        statelabel_SDV_function(filename, arrange_list);
    catch ME
        warning(['处理失败: ', filename, ', 错误信息: ', ME.message]);
    end
end

disp('批量处理完成。');


%% 
function statelabel_SDV_function(filename, arrange_list)
    bin_width = 16; % [ps]
    freq = 1.25 *10^9;  % [Hz]
    count_resol = 20;  % count resolution

    % find states of the pulses 
    period = 1 /freq *10^12;    % [ps]
    MINPEAKDISTANCE = period/bin_width - 3;

    time = csvread(filename, 0, 0);
    data = csvread(filename, 1, 0);
    [~,index_list] = findpeaks(data,'MINPEAKHEIGHT',1,'MINPEAKDISTANCE',MINPEAKDISTANCE);  
    index_list = index_list(2:end -1);
    index_first = index_list(1);
    index_last = index_list(end);

    gate_ratio = 320/800;
    pulse = zeros(1,length(index_list));
    for i = 1:length(index_list)
        pulse(i) = sum(data(index_list(i)-period /bin_width/2 * gate_ratio:1:index_list(i)+period /bin_width/2 * gate_ratio));
    end

    log_pulse = pulse;
    [y,x] = hist(log_pulse,count_resol);
    nonzero_index = find(y);
    arrset = find_continuous_sequences(nonzero_index);
    if (numel(arrset) ~= 3)
        disp(['请调整分辨率: ' filename]);
        return;
    end

    resolution = max(log_pulse)/count_resol;
    state1_range = [x(arrset{1}(1))-resolution,x(arrset{1}(end))+resolution];
    state2_range = [x(arrset{2}(1))-resolution,x(arrset{2}(end))+resolution];
    state3_range = [x(arrset{3}(1))-resolution,x(arrset{3}(end))+resolution];

    code_pulse = zeros(1,length(index_list));
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

    dic1 = struct('S',1,'D',2,'V',3);
    dic2 = struct('S',1,'D',3,'V',2);
    dic3 = struct('S',2,'D',1,'V',3);
    dic4 = struct('S',2,'D',3,'V',1);
    dic5 = struct('S',3,'D',1,'V',2);
    dic6 = struct('S',3,'D',2,'V',1);

    m = [];
    for i = 1:6
        for j = 1:length(arrange_list)
            m = [m,eval(strcat('dic',num2str(i))).(arrange_list(j))];
        end
        eval([strcat('a',num2str(i)), '= m;']);
        m = [];
    end

    for i = 1:6
        m = eval(strcat('a',num2str(i)));
        for j = 1:length(code_pulse)-length(arrange_list)
            pulse_slice = code_pulse(j:j+length(arrange_list)-1);
            pulse_slice_index = find(code_pulse(j:j+length(arrange_list)-1));
            if (m(pulse_slice_index) == pulse_slice(pulse_slice_index))
                arrange_list_order = j;
                break;
            end
        end
    end

    current_dic = eval(['dic' num2str(i)]); % 获取正确的字典
    output_data = cell(length(code_pulse)+6, 9);
    headers = {'Pulse Number', 'Pulse Counts', 'Pulse Type', 'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};

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

    pattern_counts = cell(6, 1);
    pattern_names = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};
    for pulse_idx = 1:length(code_pulse)-1
        current_type = output_data{pulse_idx, 3};
        next_type = output_data{pulse_idx+1, 3};
        next_count = output_data{pulse_idx+1, 2};
        for p = 1:6
            if strcmp(current_type, pattern_names{p}(1)) && strcmp(next_type, pattern_names{p}(2))
                pattern_counts{p} = [pattern_counts{p}; next_count];
                output_data{pulse_idx+1, p+3} = next_count;
            end
        end
    end

    stats_row = length(code_pulse)+1;
    means = zeros(1,6);
    stds = zeros(1,6);
    valid_patterns = ~cellfun(@isempty, pattern_counts);
    for p = 1:6
        if valid_patterns(p)
            means(p) = mean(pattern_counts{p});
            stds(p) = std(pattern_counts{p});
        end
    end
    ss_mean = means(1);
    sd_mean = means(4);
    norm_means = means / ss_mean;
    norm_stds = stds / ss_mean;
    ds_deviation = (means(2)-ss_mean)/ss_mean*100;
    vs_deviation = (means(3)-ss_mean)/ss_mean*100;
    dd_deviation = (means(5)-sd_mean)/sd_mean*100;
    vd_deviation = (means(6)-sd_mean)/sd_mean*100;

    output_data{stats_row, 1} = 'Mean (counts)';
    output_data{stats_row+1, 1} = 'Std Dev (counts)';
    output_data{stats_row+2, 1} = 'Normalized Mean (SS=1)';
    output_data{stats_row+3, 1} = 'Normalized Std Dev (SS=1)';
    output_data{stats_row+4, 1} = 'Deviation Percentage';

    for p = 1:6
        if valid_patterns(p)
            output_data{stats_row, p+3} = means(p);
            output_data{stats_row+1, p+3} = stds(p);
            output_data{stats_row+2, p+3} = norm_means(p);
            output_data{stats_row+3, p+3} = norm_stds(p);
        end
    end
    output_data{stats_row+4, 5} = sprintf('%.3f%% (vs SS)', ds_deviation);
    output_data{stats_row+4, 6} = sprintf('%.3f%% (vs SS)', vs_deviation);
    output_data{stats_row+4, 8} = sprintf('%.3f%% (vs SD)', dd_deviation);
    output_data{stats_row+4, 9} = sprintf('%.3f%% (vs SD)', vd_deviation);

    output_filename = strrep(filename, '.csv', '_result.xlsx');
    xlswrite(output_filename, headers, 'Sheet1', 'A1:I1');
    xlswrite(output_filename, output_data, 'Sheet1', 'A2');
end

function arrset = find_continuous_sequences(nonzero_index)
    arrset = cell(0,0);
    if isempty(nonzero_index)
        return;
    end
    start_idx = 1;
    n = numel(nonzero_index);
    while start_idx <= n
        end_idx = start_idx;
        while (end_idx < n) && (nonzero_index(end_idx)+1 == nonzero_index(end_idx+1))
            end_idx = end_idx + 1;
        end
        arrset{end+1} = nonzero_index(start_idx:end_idx);
        start_idx = end_idx + 1;
    end
end
