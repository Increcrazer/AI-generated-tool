function angle_between_polarization_states()
    % 输入第一个偏振态的斯托克斯参数 [S1, S2, S3]
    S_input1 = input('请输入第一个偏振态的斯托克斯参数 [S1, S2, S3]（例如：[1, 0, 0]）: ');
    S1 = S_input1(:);  % 确保是列向量
    
    % 输入第二个偏振态的斯托克斯参数 [S1, S2, S3]
    S_input2 = input('请输入第二个偏振态的斯托克斯参数 [S1, S2, S3]（例如：[0, 1, 0]）: ');
    S2 = S_input2(:);  % 确保是列向量
    
    % 检查输入的斯托克斯参数是否有效（归一化）
    if norm(S1) ~= 1 || norm(S2) ~= 1
        warning('输入的斯托克斯参数未归一化，程序将自动归一化');
        S1 = S1 / norm(S1);
        S2 = S2 / norm(S2);
    end
    
    % 计算两个斯托克斯向量之间的夹角（点积公式）
    dot_product = dot(S1, S2);
    angle_rad = acos(dot_product);  % 因为已经归一化，分母为1
    angle_deg = rad2deg(angle_rad);
    
    % 输出结果
    fprintf('两个偏振态对应的斯托克斯向量之间的夹角为: %.2f 度\n', angle_deg);
end