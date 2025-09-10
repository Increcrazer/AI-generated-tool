L = 0:2:100;
x_opt = zeros(length(L), 4);
fval = zeros(length(L), 1);
% raman_count = zeros(1,length(L));
raman_count = readtable('raman_count.xlsx').Total_Raman_photons_s;

if isempty(gcp('nocreate'))
    parpool('local', 40); % 使用40个工作进程
end

% 运行优化
parfor i = 1:length(L)
    [x_opt(i,:), fval(i)] = ga_qkd_optimization(L, i, raman_count);
end

% 创建表格并写入文件
headers = {'L', 'S', 'D', 'V', 'pS', 'pD', 'pV', 'SKR bit/pulse'};
data_table = [L', x_opt(:,1:2), 0.0001*ones(length(L),1), x_opt(:,3:4), 1 - x_opt(:,3) - x_opt(:,4), fval];
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
function [x_opt, fval] = ga_qkd_optimization(L, i, raman_count)
    % 参数范围设置 (k1,k2,pk1,pk2)
    nVars = 4;  % 优化变量数
    lb = [0.0001, 0.0001, 0.001, 0.001];    % 下界（避免k=0）
    ub = [1, 1, 1, 1];                % 上界
    
    % 1. 定义适应度函数
    fitnessFcn = @(x) fitnessFunction(x, L, i, raman_count);  % 负号因为GA默认最小化
    
    % 2. 定义约束条件
    Aeq = [];
    beq = [];
    
    % 非线性约束：
    %  pk3 > 0
    %  k1 > k2 + k3 
    %  k2> k3
    nonlcon = @(x) deal(...
        [ ...
         x(3) + x(4) - 1; ...
         x(2) + 0.0001 - x(1);  ...      % k1 > k2 + k3 
         0.0001 - x(2)], ...        % k2> k3
        []);                  % 无等式约束
  
    % 3. GA选项配置
    options = optimoptions('ga', ...
        'PopulationSize', 250, ...         
        'SelectionFcn', @selectionroulette, ...  
        'CrossoverFcn', @crossoverintermediate, ...   
        'CrossoverFraction', 0.8, ...      
        'MutationFcn', {@mutationadaptfeasible, 0.02}, ... 
        'MaxGenerations', 1000, ...        
        'FunctionTolerance', 1e-20, ...    
        'ConstraintTolerance', 1e-12, ...   
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

function R = fitnessFunction(x, L, i, raman_count)
    % 获取所有参数
    [epsilon_cor, epsilon_sec, qX, N, f, gate_width, width_3dB, ...
    alpha, eta_Bob_detect, e_mis_X, e_mis_Z, f_EC, ...
     p_ap, dc_count_raw, deadtime] = parameters();
    dc_count = zeros(1,4);
    for j = 1:4
        dc_count(j) = dc_count_raw(j) + raman_count(L(i)+1)/4;
    end
    k = [x(1:2) 0.0001] ;      % 前三个变量是k
    pk = [x(3:4) 1 - x(3) - x(4)];     % 后三个变量是pk
    [R_bitperpulse, ~, ~, ~, ~] = Decoy_Lim2014_corefunc(...
            epsilon_cor, epsilon_sec, ...
            qX, N, f, gate_width, width_3dB, ...
            k, pk, ...
            L(i), alpha, eta_Bob_detect, ...
            e_mis_X, e_mis_Z, f_EC, ...
            p_ap, dc_count, deadtime);
    R = -R_bitperpulse; % 取负因为GA默认最小化
end

function [epsilon_cor, epsilon_sec, qX, N, f, gate_width, width_3dB, ...
         alpha, eta_Bob_detect, e_mis_X, e_mis_Z, f_EC, ...
          p_ap, dc_count_raw, deadtime] = parameters()
    % 所有参数定义
    epsilon_cor = 10^-15;
    epsilon_sec = 10^-15;
    N = 10^13;
    qX = 0.5;
    f = 1.25*10^9;
    gate_width = 320*10^(-12);
    width_3dB = 70*10^(-12);
    alpha = 0.2;
    eta_Bob_detect = 0.2;
    e_mis_X = 0.01;
    e_mis_Z = 0.01;
    f_EC = 1.16;
    p_ap = 0.04;
    deadtime = 50*10^(-9);
    dc_count_raw = [125 125 125 125];
end
