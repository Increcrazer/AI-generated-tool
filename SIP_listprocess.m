% 将目录下所有xx_raw.txt转化成xx.txt文件，输入给prbs_generator文件
% 获取当前目录下所有xx_raw.txt文件
files = dir('*_raw.txt');

% 对每个文件进行处理
for k = 1:length(files)
    % 读取原始文件名
    raw_filename = files(k).name;
    
    % 确定输出文件名（去掉_raw）
    output_filename = strrep(raw_filename, '_raw', '');
    
    % 读取文件内容 - 使用更可靠的方式
    str = fileread(raw_filename);
    str = str(1:256);  % 确保只取前256个字符
    
    % 检查字符数并显示警告
    if length(str) > 256
        fprintf('警告: 文件 %s 包含 %d 个字符，将只使用前256个字符\n', raw_filename, length(str));
    elseif length(str) < 256
        error('文件 %s 不足256个字符，实际有%d个字符', raw_filename, length(str));
    end
    
    % 1. 将序列分成32行8列的矩阵
    matrix_32x8 = reshape(str, 8, 32)';
    
    % 2. 转置矩阵为8行32列
    matrix_8x32 = matrix_32x8';
    
    % 3. 列对调：第一列和最后一列对调，第二列和倒数第二列对调，以此类推
    cols = size(matrix_8x32, 2);
    for i = 1:floor(cols/2)
        temp = matrix_8x32(:, i);
        matrix_8x32(:, i) = matrix_8x32(:, cols - i + 1);
        matrix_8x32(:, cols - i + 1) = temp;
    end
    
    % 4. 在每个字符后面加上'F'
    final_matrix = repmat(' ', 8, 32*2); % 预分配空间
    for row = 1:8
        for col = 1:32
            final_matrix(row, (col-1)*2+1) = matrix_8x32(row, col);
            final_matrix(row, (col-1)*2+2) = 'F';
        end
    end
    
    % 将结果写入新文件
    fid = fopen(output_filename, 'w');
    for row = 1:8
        fprintf(fid, '%s\n', final_matrix(row, :));
    end
    fclose(fid);
    
    fprintf('已处理文件: %s → %s\n', raw_filename, output_filename);
end

disp('所有文件处理完成！');