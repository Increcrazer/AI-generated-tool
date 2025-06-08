%% user set
filename = 'SDV-IM.csv';
bin_width = 16;
bin_number = 100000;
freq = 1.25 *10^9;  % [Hz]
pseudo_length = 32;
count_resol = 30;  % count resolution
arrange_list = char("SDVVSSDDSSVSVSSDSSVD");

%% find states of the pulses 
period = 1 /freq *10^12;    % [ps]
pseudo_time = 3*pseudo_length *period;   % [ps] 取3倍伪随机数长度

time = csvread(filename, 0, 0,[0,0,0,pseudo_time /bin_width]);
data = csvread(filename, 1, 0,[1,0,1,pseudo_time /bin_width]);
[~,index_first] = max(data(1:5*period /bin_width)); % 取5个脉冲内的第一个峰值
index_last = index_first + (pseudo_time -5*period) /bin_width;
index_list = index_first:period /bin_width:index_last;

gate_ratio = 320/800;
pulse = zeros(1,length(index_list));
for i = 1:length(index_list)
    pulse(i) = sum(data(index_list(i)-period /bin_width/2 * gate_ratio:1:index_list(i)+period /bin_width/2 * gate_ratio));
end

log_pulse = log10(pulse);
figure(1);
[y,x] = hist(log_pulse,count_resol);   % y：每个柱子中包含的数据点数量（频数） x：每个柱子的中心位置坐标
bar(x,y);
xlabel('counts/log','FontName','Times New Roman','fontsize',18);
ylabel('probability','FontName','Times New Roman','fontsize',18);
title('Histogram');
nonzero_index = find(y); % find 函数用于查找数组中非零元素的索引

arrset = find_continuous_sequences(nonzero_index);

if (numel(arrset) ~= 3)
    disp('请调整分辨率');
end

resolution = max(log_pulse)/count_resol;
state1_range = [x(arrset{1}(1))-resolution,x(arrset{1}(end))+resolution];
state2_range = [x(arrset{2}(1))-resolution,x(arrset{2}(end))+resolution];
state3_range = [x(arrset{3}(1))-resolution,x(arrset{3}(end))+resolution];

% code the pulses
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

% 实际上，SDV的情况毋庸置疑是dic6，不是说明有问题
% find pulse position
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
            disp('Find position!');
            disp(strcat('Dictionary is dic',num2str(i)));
            disp(strcat('a(1) corresponds to code_pulse(',num2str(j),')'));
            arrange_list_order = j;
            break;
        else
            continue;
        end
    end
end

%% 将pulse序列归类并写入Excel，计算完整统计量并格式化输出
% 使用上一步找到的dic{i}和起始位置arrange_list_order

% 创建输出矩阵（增加统计行）
output_data = cell(length(code_pulse)+6, 9); % 增加6行用于统计结果

% 列标题
headers = {'Pulse Number', 'Pulse Counts', 'Pulse Type', 'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};

% 使用找到的dic{i}进行归类
current_dic = eval(['dic' num2str(i)]); % 获取正确的字典

% 填充前三列数据
for pulse_idx = 1:length(code_pulse)
    % 第一列：脉冲序号
    output_data{pulse_idx, 1} = pulse_idx;
    
    % 第二列：脉冲计数值
    output_data{pulse_idx, 2} = pulse(pulse_idx);
    
    % 第三列：脉冲类型
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

% 初始化模式检测标志和数据收集
pattern_detected = false(length(code_pulse), 6);
pattern_counts = cell(6, 1);
pattern_names = {'SS', 'DS', 'VS', 'SD', 'DD', 'VD'};

% 检测模式并收集数据
for pulse_idx = 1:length(code_pulse)-1
    current_type = output_data{pulse_idx, 3};
    next_type = output_data{pulse_idx+1, 3};
    next_count = output_data{pulse_idx+1, 2};
    
    for p = 1:6
        first_char = pattern_names{p}(1);
        second_char = pattern_names{p}(2);
        
        if strcmp(current_type, first_char) && strcmp(next_type, second_char)
            pattern_detected(pulse_idx+1, p) = true;
            pattern_counts{p} = [pattern_counts{p}; next_count];
            output_data{pulse_idx+1, p+3} = next_count;
        end
    end
end

% 计算统计量 ====================================================
stats_row = length(code_pulse)+1;

% 计算基本统计量
means = zeros(1,6);
stds = zeros(1,6);
valid_patterns = ~cellfun(@isempty, pattern_counts);

for p = 1:6
    if valid_patterns(p)
        means(p) = mean(pattern_counts{p});
        stds(p) = std(pattern_counts{p});
    end
end

% 归一化因子
ss_mean = means(1);
sd_mean = means(4);

% 计算归一化值和偏差百分比
norm_means = means / ss_mean;
norm_stds = stds / ss_mean;

% DS和VS相对于SS的偏差百分比
if valid_patterns(2)
    ds_deviation = (means(2)-ss_mean)/ss_mean*100;
else
    ds_deviation = NaN;
end

if valid_patterns(3)
    vs_deviation = (means(3)-ss_mean)/ss_mean*100;
else
    vs_deviation = NaN;
end

% DD和VD相对于SD的偏差百分比
if valid_patterns(5) && valid_patterns(4)
    dd_deviation = (means(5)-sd_mean)/sd_mean*100;
else
    dd_deviation = NaN;
end

if valid_patterns(6) && valid_patterns(4)
    vd_deviation = (means(6)-sd_mean)/sd_mean*100;
else
    vd_deviation = NaN;
end

% 填充统计行 ===================================================
% 平均值行
output_data{stats_row, 1} = 'Mean (counts)';
for p = 1:6
    if valid_patterns(p)
        output_data{stats_row, p+3} = means(p);
    end
end

% 标准差行
output_data{stats_row+1, 1} = 'Std Dev (counts)';
for p = 1:6
    if valid_patterns(p)
        output_data{stats_row+1, p+3} = stds(p);
    end
end

% 归一化均值行
output_data{stats_row+2, 1} = 'Normalized Mean (SS=1)';
for p = 1:6
    if valid_patterns(p)
        output_data{stats_row+2, p+3} = norm_means(p);
    end
end

% 归一化标准差行
output_data{stats_row+3, 1} = 'Normalized Std Dev (SS=1)';
for p = 1:6
    if valid_patterns(p)
        output_data{stats_row+3, p+3} = norm_stds(p);
    end
end

% 偏差百分比行
output_data{stats_row+4, 1} = 'Deviation Percentage';
if ~isnan(ds_deviation)
    output_data{stats_row+4, 5} = sprintf('%.3f%% (vs SS)', ds_deviation);
end
if ~isnan(vs_deviation)
    output_data{stats_row+4, 6} = sprintf('%.3f%% (vs SS)', vs_deviation);
end
if ~isnan(dd_deviation)
    output_data{stats_row+4, 8} = sprintf('%.3f%% (vs SD)', dd_deviation);
end
if ~isnan(vd_deviation)
    output_data{stats_row+4, 9} = sprintf('%.3f%% (vs SD)', vd_deviation);
end

% 写入Excel文件并设置格式 ======================================
output_filename = strrep(filename, '.csv', '_result.xlsx')

% 先写入数据
xlswrite(output_filename, headers, 'Sheet1', 'A1:I1');
xlswrite(output_filename, output_data, 'Sheet1', 'A2');

% 设置Excel格式（需要Excel COM接口）
try
    excel = actxserver('Excel.Application');
    workbook = excel.Workbooks.Open(fullfile(pwd, output_filename));
    sheet = workbook.Sheets.Item(1);
    sheet.Columns.AutoFit;
    
    % 设置所有数据居中
    all_range = sheet.Range('A1:I'+string(size(output_data,1)+1));
    all_range.HorizontalAlignment = 3; % xlCenter
    
    % 设置统计行加粗
    for row = stats_row:stats_row+4
        sheet.Range('A'+string(row+1)+':I'+string(row+1)).Font.Bold = true;
        num_range = sheet.Range('A'+string(row+1)+':I'+string(row+1));
        num_range.NumberFormat = '0.0000';
    end
    
    % 保存并关闭
    workbook.Save;
    workbook.Close;
    excel.Quit;
    
catch
    warning('自动格式设置失败，请手动调整Excel格式');
end

disp('===== 关键统计结果 =====');
disp(['SS均值: ' num2str(ss_mean) ' (归一化基准)']);
disp(['SD均值: ' num2str(sd_mean)]);
disp(['DS相对于SS偏差: ' num2str(ds_deviation, '%.3f') '%']);
disp(['VS相对于SS偏差: ' num2str(vs_deviation, '%.3f') '%']);
disp(['DD相对于SD偏差: ' num2str(dd_deviation, '%.3f') '%']);
disp(['VD相对于SD偏差: ' num2str(vd_deviation, '%.3f') '%']);
disp(['结果已写入文件: ' output_filename]);

%%  plot and label state    
figure(2);
time = time(index_first-period /bin_width /2:index_last-period /bin_width /2);
data = data(index_first-period /bin_width /2:index_last-period /bin_width /2);
H = semilogy(time,data);
grid on;
xmin = time((arrange_list_order-1) *period /bin_width);
xmax = time((arrange_list_order+length(arrange_list)-1) *period /bin_width);
xlim([xmin, xmax]);
% ylim([0,5]);
xlabel('time/ps','FontName','Times New Roman','fontsize',18);
ylabel('counts/log','FontName','Times New Roman','fontsize',18);
title('SDV');
xpos = xmin+period/2:period:xmax-period/2;
ypos = 8000;
for i = 1:length(xpos)
    text(xpos(i),ypos,arrange_list(i));
end

%% 连续数字序列查找函数
function arrset = find_continuous_sequences(nonzero_index)
    arrset = cell(0,0);
    if isempty(nonzero_index)
        return;
    end
    
    start_idx = 1;
    n = numel(nonzero_index);
    
    while start_idx <= n
        % 查找当前连续序列的结束位置
        end_idx = start_idx;
        while (end_idx < n) && (nonzero_index(end_idx)+1 == nonzero_index(end_idx+1))
            end_idx = end_idx + 1;
        end
        
        % 保存当前连续序列
        arrset{end+1} = nonzero_index(start_idx:end_idx);
        
        % 移动到下一个序列的起始位置
        start_idx = end_idx + 1;
    end
end
