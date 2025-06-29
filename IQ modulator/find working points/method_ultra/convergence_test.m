% 运行1000次并统计迭代次数分布
num_runs = 1000;
iteration_results = zeros(num_runs, 1);

for i = 1:num_runs
    iteration_results(i) = find_IQ_phase_method_ultra_with_count();
end

% 显示统计结果
fprintf('迭代次数统计（共%d次运行）：\n', num_runs);
fprintf('平均迭代次数: %.2f\n', mean(iteration_results));
fprintf('最小迭代次数: %d\n', min(iteration_results));
fprintf('最大迭代次数: %d\n', max(iteration_results));

% 绘制直方图
figure;
histogram(iteration_results, 'BinWidth', 1);
title('迭代次数分布');
xlabel('迭代次数');
ylabel('出现频率');
grid on;
%% 计算迭代次数
function iteration_counts = find_IQ_phase_method_ultra_with_count()
    % 初始参数设置
    IQ_IN = 1;
    R = 0.5;
    I_RF_phi = 0;
    Q_RF_phi = 0;
    k = -1.5;
    b = 7;

    % 初始随机相位 (0-2pi)
    IQ_fai0 = 2*pi*rand();
    up_fai0 = 2*pi*rand();
    down_fai0 = 2*pi*rand();
    
    % 初始相位偏移
    PS2_fai = 0;
    PS3_fai = 0;
    PS5_fai = 0;

    % 最大迭代次数
    max_iterations_1 = 20;
    max_iterations_23 = 1;

    % 初始化前一次PS5功率和PS2功率
    prev_PS5_power = -Inf;
    
    % 全局循环终止标志
    global_converged = false;
    
    % 迭代计数器
    total_iterations = 0;
        
    for iter_1 = 1:max_iterations_1  
        if global_converged
            break;
        end
        
        % 扫描PS5_fai
        [PS5_fai, PS5_power] = scan_phase_mid(PS2_fai, PS3_fai, PS5_fai, ...
                                         IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                         I_RF_phi, Q_RF_phi, k, b, 'PS5');
        
        % 从第二次迭代开始比较功率
        if iter_1 > 1 && PS5_power > prev_PS5_power
            PS5_fai = mod(prev_PS5_fai + pi, 2*pi);  % 加π并取模
        end
        prev_PS5_fai = PS5_fai;
        prev_PS5_power = PS5_power;
        
        for iter_23 = 1:max_iterations_23    
            if global_converged
                break;
            end
            % 扫描PS2_fai
            [PS2_fai, ~] = scan_phase_min(PS2_fai, PS3_fai, PS5_fai, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, 'PS2');
            PS2_fai = mod(PS2_fai,2*pi);

            % 扫描PS3_fai
            [PS3_fai, ~] = scan_phase_min(PS2_fai, PS3_fai, PS5_fai, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, 'PS3');
            PS3_fai = mod(PS3_fai,2*pi);
            
            % 计算组合相位
            combined_PS2 = mod(PS2_fai + up_fai0, 2*pi);
            combined_PS3 = mod(PS3_fai + down_fai0, 2*pi);
            
            % 检查是否接近π
            if abs(combined_PS2 - pi)/pi<1/300 && abs(combined_PS3 - pi)/pi<1/300
                global_converged = true;
            end 
        end
        % 更新迭代计数器
        total_iterations = total_iterations + 1;
    end
    % 返回总迭代次数
    iteration_counts = total_iterations;
end

%% 辅助函数：相位扫描
function [best_phase, mid_power] = scan_phase_mid(PS2, PS3, PS5, ...
                                         IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                         I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 600); % 扫描600个点
    power_values = zeros(size(phase_points));
    
    % 首先收集所有相位点的功率值
    for i = 1:length(phase_points)
        phase = phase_points(i);
        
        % 根据目标更新相应的相位
        switch target
            case 'PS2'
                current_PS2 = phase;
                current_PS3 = PS3;
                current_PS5 = PS5;
            case 'PS3'
                current_PS2 = PS2;
                current_PS3 = phase;
                current_PS5 = PS5;
            case 'PS5'
                current_PS2 = PS2;
                current_PS3 = PS3;
                current_PS5 = phase;
        end
        
        % 计算光功率
        I = abs(layout_IQ(IQ_IN, R, IQ_fai0, ...
             up_fai0, I_RF_phi, current_PS2, ...
             down_fai0, Q_RF_phi, current_PS3, ...
             current_PS5, k, b))^2;
        power_values(i) = I;
    end
    
    % 计算最大和最小功率
    max_power = max(power_values);
    min_power = min(power_values);
    target_power = (max_power + min_power)/2;
    
    % 找到最接近目标功率的相位点
    [~, idx] = min(abs(power_values - target_power));
    best_phase = phase_points(idx);
    mid_power = power_values(idx);
end

function [best_phase, min_power] = scan_phase_min(PS2, PS3, PS5, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 600); % 扫描600个点
    min_power = inf;
    best_phase = eval(target); % 获取当前目标相位值
    
    for phase = phase_points
        % 根据目标更新相应的相位
        switch target
            case 'PS2'
                current_PS2 = phase;
                current_PS3 = PS3;
                current_PS5 = PS5;
            case 'PS3'
                current_PS2 = PS2;
                current_PS3 = phase;
                current_PS5 = PS5;
            case 'PS5'
                current_PS2 = PS2;
                current_PS3 = PS3;
                current_PS5 = phase;
        end
        
        % 计算光功率
        I = abs(layout_IQ(IQ_IN, R, IQ_fai0, ...
             up_fai0, I_RF_phi, current_PS2, ...
             down_fai0, Q_RF_phi, current_PS3, ...
             current_PS5, k, b))^2;
        current_power = I;
        
        if current_power < min_power
            min_power = current_power;
            best_phase = phase;
        end
    end
end

%% 2to2 IQ in layout
function Eo = layout_IQ(IQ_IN, R, IQ_fai0, ...
    up_fai0, I_RF_phi, PS2_fai, ...
    down_fai0, Q_RF_phi, PS3_fai, ...
    PS5_fai, k, b)
    [~,Eo] = IQ(IQ_IN,0,R, IQ_fai0, ...
    up_fai0,0,I_RF_phi,PS2_fai,0, ... 
    down_fai0,0,Q_RF_phi,PS3_fai,0, ... 
    PS5_fai,0, ...
    k,b);
% down_fai0和PS3_fai加在同意臂，是为了直观显示IQ自动搜索消光点算法
end

%% 2to2 IQ 
function [Eo1,Eo2] = IQ (Ein1,Ein2,R, IQ_fai0, ...
up_fai0,up_CDM_fai1,up_CDM_fai2,up_PS_fai1,up_PS_fai2, ... 
down_fai0,down_CDM_fai1,down_CDM_fai2,down_PS_fai1,down_PS_fai2, ... 
IQ_up_PS_fai,IQ_down_PS_fai, ...
k,b)
    [E1,E2] = BS(R,Ein1,Ein2);
    [~,E3] = MZI(E1,0,R,up_fai0,up_CDM_fai1,up_CDM_fai2,up_PS_fai1,up_PS_fai2,k,b);
    [~,E4] = MZI(E2,0,R,down_fai0,down_CDM_fai1,down_CDM_fai2,down_PS_fai1,down_PS_fai2,k,b);
    E5 = E3*exp(1i*IQ_up_PS_fai)*exp(1i*IQ_fai0);
    E6 = E4*exp(1i*IQ_down_PS_fai);
    [Eo1,Eo2] = BS(R,E5,E6);
end

%% 2to2 MZI 
function [Eo1,Eo2] = MZI(Ein1,Ein2,R,fai0,CDM_fai1,CDM_fai2,PS_fai1,PS_fai2,k,b)
    % fai_0是上下臂固有相位差
    loss_factor1 = k*CDM_fai1+b;   % k is negative
    loss_factor2 = k*CDM_fai2+b;
    [E1,E2] = BS(R,Ein1,Ein2);
    E3 = E1*exp(1i*CDM_fai1)*sqrt(10^(-loss_factor1/10));
    E4 = E2*exp(1i*CDM_fai2)*sqrt(10^(-loss_factor2/10));
    E5 = E3*exp(1i*PS_fai1)*exp(1i*fai0);
    E6 = E4*exp(1i*PS_fai2);
    [Eo1,Eo2] = BS(R,E5,E6);
end

%% 2to2 BS
function [Eo1,Eo2] = BS(R,Ein1,Ein2)
    T = 1-R;
    Eo1 = 1i*sqrt(R)*Ein1 + sqrt(T)*Ein2;
    Eo2 = sqrt(T)*Ein1 + 1i*sqrt(R)*Ein2;
end