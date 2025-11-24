% 生成256个随机字符（'P', 'R', 'N', 'L'），概率比为1:1:1:1
chars = ['P', 'R', 'N', 'L'];
probabilities = [1, 1, 1, 1]; % P:R:N:L = 1:1:1:1
cum_prob = cumsum(probabilities) / sum(probabilities);
random_seq = chars(arrayfun(@(x) find(x <= cum_prob, 1), rand(1, 1024)));

% 保存到randsource.txt
fid = fopen('randsource_1024.txt', 'w');
fprintf(fid, '%s', random_seq);
fclose(fid);

% 定义映射规则
mapping.PM1 = struct('P', '00', 'N', 'F0', 'R', '00', 'L', 'F0');
mapping.PM2 = struct('P', '00', 'N', 'F0', 'R', 'F0', 'L', 'F0');
mapping.PM3 = struct('P', '00', 'N', '00', 'R', '00', 'L', 'F0');

% 生成PM1.txt, PM2.txt
pm_files = {'PM1_1024.txt', 'PM2_1024.txt', 'PM3_1024.txt'};
pm_names = {'PM1_1024', 'PM2_1024', 'PM3_1024'};

for i = 1:3
    pm_name = pm_names{i};
    mapped_seq = '';
    
    % 根据映射规则转换字符
    for c = random_seq
        switch c
            case 'P'
                mapped_seq = [mapped_seq, mapping.(pm_name).P];
            case 'N'
                mapped_seq = [mapped_seq, mapping.(pm_name).N];
            case 'R'
                mapped_seq = [mapped_seq, mapping.(pm_name).R];
            case 'L'
                mapped_seq = [mapped_seq, mapping.(pm_name).L];
        end
    end
    
    % 写入文件（带BOM头）
    fid = fopen(pm_files{i}, 'w');
    fwrite(fid, [239 187 191], 'uint8'); % UTF-8 BOM
    fprintf(fid, '%s', mapped_seq);
    fclose(fid);
end

disp('文件生成完成：randsource_1024.txt, PM1_1024.txt, PM2_1024.txt, PM3_1024.txt');
