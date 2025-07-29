% 生成1024个随机字符（'S', 'D', 'V'），概率比为2:1:1
chars = ['S', 'D', 'V'];
probabilities = [2, 1, 1]; % S:D:V = 2:1:1
cum_prob = cumsum(probabilities) / sum(probabilities);
random_seq = chars(arrayfun(@(x) find(x <= cum_prob, 1), rand(1, 1024)));

% 保存到randsource.txt
fid = fopen('randsource.txt', 'w');
fprintf(fid, '%s', random_seq);
fclose(fid);

% 定义映射规则
mapping.AM1 = struct('S', 'F0', 'D', 'F0', 'V', '00');
mapping.AM2 = struct('S', 'F0', 'D', '00', 'V', '00');

% 生成AM1.txt, AM2.txt
pm_files = {'AM1.txt', 'AM2.txt'};
pm_names = {'AM1', 'AM2'};

for i = 1:2
    pm_name = pm_names{i};
    mapped_seq = '';
    
    % 根据映射规则转换字符
    for c = random_seq
        switch c
            case 'S'
                mapped_seq = [mapped_seq, mapping.(pm_name).S];
            case 'D'
                mapped_seq = [mapped_seq, mapping.(pm_name).D];
            case 'V'
                mapped_seq = [mapped_seq, mapping.(pm_name).V];
        end
    end
    
    % 写入文件（带BOM头）
    fid = fopen(pm_files{i}, 'w');
    fwrite(fid, [239 187 191], 'uint8'); % UTF-8 BOM
    fprintf(fid, '%s', mapped_seq);
    fclose(fid);
end

disp('文件生成完成：randsource.txt, AM1.txt, AM2.txt');
