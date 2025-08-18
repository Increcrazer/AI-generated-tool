L = 0:10:250;
x_opt = zeros(length(L), 6);
fval = zeros(length(L), 1);

% 运行优化
for i = 1:length(L)
    [x_opt(i,:), fval(i)] = ga_qkd_optimization(L(i));
end

% 创建表格并写入文件
headers = {'L', 'S', 'D', 'V', 'pS', 'pD', 'pV', 'SKR bit/pulse'};
data_table = [L', x_opt, fval];
T = array2table(data_table, 'VariableNames', headers);

% 保存为CSV文件
writetable(T, 'optimization_results.csv');

% Create figure with improved styling
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% Plot with logarithmic y-axis
semilogy(L, fval, 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]);
grid on;

% Add labels and title with larger fonts
xlabel('Distance (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Secret Key Rate (bits/pulse)', 'FontSize', 14, 'FontWeight', 'bold');
title('QKD Key Rate vs Distance', 'FontSize', 16, 'FontWeight', 'bold');

% Adjust axes properties
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.YScale = 'log';
ax.YMinorTick = 'on';
ax.XMinorTick = 'on';
ax.GridAlpha = 0.3;

%%
function [x_opt, fval] = ga_qkd_optimization(L)
    % 参数范围设置 (k1,k2,k3,pk1,pk2,pk3)
    nVars = 6;  % 优化变量数
    lb = [0.0001, 0.0001, 0.0001, 0.001, 0.001, 0.001];    % 下界（避免k=0）
    ub = [1, 1, 1, 1, 1, 1];                % 上界
    
    % 1. 定义适应度函数
    fitnessFcn = @(x) fitnessFunction(x, L);  % 负号因为GA默认最小化
    
    % 2. 定义约束条件
    % 线性等式约束：pk1 + pk2 + pk3 = 1
    Aeq = [0, 0, 0, 1, 1, 1];
    beq = 1;
    
    % 非线性约束：
    % (1) k1/k3 > 100 
    % (2) k1/k2 < 10
    % (3) k1/k2 >1.25
    % (3) k1 > k2 > k3
    nonlcon = @(x) deal(...
        [-(x(1)/x(3) - 100); ...  % k1/k3 > 100 → c1 ≤ 0
         (x(1)/x(2) - 10); ... % k1/k2 < 10 → c2 ≤ 0
         -(x(1)/x(2) - 1.25); ... % k1/k2 >1.25 → c3 ≤ 0
         x(2) - x(1);  ...      % k1 > k2 → c4 ≤ 0
         x(3) - x(2)], ...        % k2 > k3 → c5 ≤ 0
        []);                  % 无等式约束
  
    % 3. GA选项配置
    options = optimoptions('ga', ...
        'PopulationSize', 150, ...         
        'SelectionFcn', @selectionroulette, ...  
        'CrossoverFcn', @crossovertwopoint, ...   
        'CrossoverFraction', 0.8, ...      
        'MutationFcn', {@mutationadaptfeasible, 0.02}, ... 
        'MaxGenerations', 1000, ...        
        'FunctionTolerance', 1e-20, ...    
        'ConstraintTolerance', 1e-10, ...   
        'StallGenLimit', 50, ...       
        'Display', 'iter', ...             
        'PlotFcn', {@gaplotbestf, @gaplotstopping} ...
        );
    
    % 4. 运行GA优化（带约束）
    [x_opt, fval, ~] = ga(fitnessFcn, nVars, ...
        [], [], Aeq, beq, lb, ub, nonlcon, options);
    fval = -fval;
    
%     % 5. 结果输出
%     optimal_params = struct(...
%         'k', x_opt(1:3), ...
%         'pk', x_opt(4:6), ...
%         'key_rate', fval, ...
%         'constraint_violation', max([...
%             0, ...
%             x_opt(1)-x_opt(2), ...       % k1 > k2
%             x_opt(2)-x_opt(3) ...        % k2 > k3
%         ]));
%     
%     disp('优化结果:');
%     disp(['k = [', num2str(optimal_params.k, '%.4f '), ']']);
%     disp(['pk = [', num2str(optimal_params.pk, '%.4f '), ']']);
%     disp(['密钥率 = ', num2str(optimal_params.key_rate)]);
%     disp(['最大约束违反量 = ', num2str(optimal_params.constraint_violation)]);
end

function R = fitnessFunction(x, L)
    % 获取所有参数
    [epsilon_cor, epsilon_sec, qX, N, f, gate_width, width_3dB, ...
    alpha, eta_Bob_detect, e_mis_X, e_mis_Z, f_EC, ...
     p_ap, dc_count, deadtime] = parameters();
    
    k = x(1:3);      % 前三个变量是k
    pk = x(4:6);     % 后三个变量是pk
    [R_bitperpulse, ~, ~, ~, ~] = Decoy_Lim2014_corefunc(...
            epsilon_cor, epsilon_sec, ...
            qX, N, f, gate_width, width_3dB, ...
            k, pk, ...
            L, alpha, eta_Bob_detect, ...
            e_mis_X, e_mis_Z, f_EC, ...
            p_ap, dc_count, deadtime);
    R = -R_bitperpulse; % 取负因为GA默认最小化
end

function [c, ceq] = constraints(x)
    % 提取变量
    k = x(1:3);  % k(1), k(2), k(3)
    pk = x(4:6); % pk(1), pk(2), pk(3)
    
    % 非线性不等式约束：k₁/k₃ > 100 → 转换为 k₁ - 100*k₃ > 0
    c = -(k(1) - 100*k(3)); % GA默认处理 c ≤ 0，因此取负
    
    % 非线性等式约束（无，保留原线性等式约束）
    ceq = [];
end

function [epsilon_cor, epsilon_sec, qX, N, f, gate_width, width_3dB, ...
         alpha, eta_Bob_detect, e_mis_X, e_mis_Z, f_EC, ...
          p_ap, dc_count, deadtime] = parameters()
    % 所有参数定义
    epsilon_cor = 10^-15;
    epsilon_sec = 10^-15;
    N = 10^13;
    qX = 0.5;
    f = 1.25*10^9;
    gate_width = 1/f;
    width_3dB = 70*10^(-12);
    alpha = 0.2;
    eta_Bob_detect = 0.4;
    e_mis_X = 0.005;
    e_mis_Z = 0.005;
    f_EC = 1.16;
    p_ap = 0.04;
    dc_count = [100 100 100 100];
    deadtime = 0;
end
