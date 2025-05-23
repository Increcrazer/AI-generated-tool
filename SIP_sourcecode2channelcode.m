% 将xx.txt内的码字流按照通道名称进行转换，方便示波器上进行核对。
% 转换映射关系由'编码关系.xlsx'给出
% 默认有8个GTX通道

% 1. 读取编码关系表（确保码字为文本格式，且无空格）
[~, ~, code_table] = xlsread('编码关系.xlsx');
code_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 2:size(code_table, 1)
    % 处理码字：转为大写、去除空格、强制文本格式
    hex_code = upper(strtrim(num2str(code_table{i, 1})));
    outputs = cell2mat(code_table(i, 2:9));
    code_map(hex_code) = outputs;
end

% 2. 处理_raw.txt文件
raw_files = dir('*_raw.txt');
for file_idx = 1:length(raw_files)
    raw_filename = raw_files(file_idx).name;
    file_prefix = raw_filename(1:strfind(raw_filename, '_raw')-1);
    
    % 读取文件并提取1位十六进制码字
    file_content = fileread(raw_filename);
    file_content = file_content(1:end-1);
    hex_codes = regexp(file_content, '[0-9A-Fa-f]', 'match');
    
    % 初始化通道数据
    channels = { 'XX','PM2', 'PM1', 'PM0', 'AM2', 'AM1', 'AM0', 'LD'};
    channel_data = struct();
    for ch = channels
        channel_data.(ch{1}) = [];
    end
    
    % 3. 译码
    for i = 1:length(hex_codes)
        hex_code = upper(hex_codes{i});
        if isKey(code_map, hex_code)
            outputs = code_map(hex_code);
            for ch_idx = 1:length(channels)
                 channel_data.(channels{ch_idx}) = [channel_data.(channels{ch_idx}), outputs(ch_idx)];
            end
        else
            warning('码字 %s 在编码表中未找到，请检查Excel文件！', hex_code);
        end
    end
    
    % 4. 保存文件
    for ch = channels
        ch_name = ch{1};
        data = channel_data.(ch_name);
        if ~isempty(data)
            output_filename = sprintf('%s_%s.txt', file_prefix, ch_name);
            fid = fopen(output_filename, 'w');
            fprintf(fid, '%d', data);
            fclose(fid);
        end
    end
end
disp('处理完成。');