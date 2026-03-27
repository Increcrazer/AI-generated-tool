%% ===================== 1. 定义输入数据 =====================
hex_strings = {
    'LNLLPLNRRRLRNPLLPRNNLLPNPRPNRNRN';
    'RNLRRPRRPNRLRPLPNLPLRPPRPLPLLNLL';
    'PRPRNRNPRLNRLRPNNPLRRPNLLNRRPPNP';
    'PLLNRPPNPPRNRNLRNPRLPLPLNNRRRPRN';
    'PLLRPNRPNRNLNNNRPLPPRLPPRPLNLNLL';
    'NNPLRNPRLLNRNRNLNLPPNPLLRNNPNLLL';
    'LLRPPRLLLNLRNRPNNLNPRLNLNNPNNPNL';
    'PLPPRLPRRRRLLPLPPPLPRPPLLLPNLRLL'
};

%% ===================== 2. 先做倒序+取偶数列处理 =====================
result = cell(size(hex_strings));
for i = 1:length(hex_strings)
    current_str = char(hex_strings{i}); 
    reversed_str = fliplr(current_str);
    final_str = reversed_str(1:1:end); 
    result{i} = final_str;
end

%% ===================== 3. 按列提取合成大字符串 =====================
num_rows = length(result);
num_cols = length(result{1});
for i = 2:num_rows
    if length(result{i}) ~= num_cols
        error('错误：第%d行长度与第1行不一致，无法按列提取！', i);
    end
end

big_string = ''; 
for col = 1:num_cols          
    for row = 1:num_rows      
        big_string = [big_string, result{row}(col)];
    end
end

%% ===================== 4. 【核心新增】字符映射替换 =====================
% --- 4.1 统一转为大写，避免小写a-f匹配不上 ---
big_string_upper = upper(big_string);


%% ===================== 5. 输出结果 =====================
fprintf('========== 处理完成 ==========\n');
fprintf('按列合成的大字符串长度：%d\n', length(big_string));

% 输出前100个字符对比
fprintf('--- 前256个字符对比 ---\n');
fprintf('原大字符串: %s\n', big_string(1:min(256, end)));
