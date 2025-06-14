function find_IQ_phase_method_ultra()
    % 一种找到I和Q消光点的迭代算法（无敌）
    % 主函数封装所有变量
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
    
    fprintf('up_fai0 = %.4f rad, down_fai0 = %.4f rad, IQ_fai0 = %.4f rad\n', ...
    up_fai0, down_fai0, IQ_fai0);

    % 初始相位偏移
    PS2_fai = 0;
    PS3_fai = 0;
    PS5_fai = 0;

    % 迭代次数
    num_iterations_1 = 10;
    num_iterations_23 = 1;

    % 存储结果
    results = [];
    result_count = 0;
    
    % 初始化前一次PS5功率和PS2功率
    prev_PS5_power = -Inf;
    prev_PS2_power = -Inf;
    
    % 添加全局循环终止标志
    global_converged = false;
        
    for iter_1 = 1:num_iterations_1  
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
        % 更新前一次功率记录
        prev_PS5_power = PS5_power;
        
        for iter_23 = 1:num_iterations_23    
            if global_converged
                break;
            end
            
            % 保存旧的功率值用于比较
            old_PS2_power = prev_PS2_power;

            % 扫描PS2_fai
            [PS2_fai, PS2_power] = scan_phase_min(PS2_fai, PS3_fai, PS5_fai, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, 'PS2');
            PS2_fai = mod(PS2_fai,2*pi);

            % 扫描PS3_fai
            [PS3_fai, PS3_power] = scan_phase_min(PS2_fai, PS3_fai, PS5_fai, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, 'PS3');
            PS3_fai = mod(PS3_fai,2*pi);
            
            % 计算组合相位
            combined_PS2 = mod(PS2_fai + up_fai0, 2*pi);
            combined_PS3 = mod(PS3_fai + down_fai0, 2*pi);
            
            % 检查是否接近π
            if abs(combined_PS2 - pi) < 0.01 && abs(combined_PS3 - pi) < 0.01
                fprintf('全局收敛条件满足：PS2_fai+up_fai0=%.4f, PS3_fai+down_fai0=%.4f 都接近π\n',...
                        combined_PS2, combined_PS3);
                global_converged = true;
            end
            
            % 记录结果
            result_count = result_count + 1;
            results(result_count,:) = [iter_1, iter_23, PS2_fai, PS3_fai, PS5_fai, ...
                                      PS2_power, PS3_power, PS5_power];
            
            % 检查PS2功率是否收敛（变化小于0.0001）
            if iter_23 > 1 && abs(PS2_power - old_PS2_power)/old_PS2_power < 0.001
                break;
            end

            % 更新前一次PS2功率记录
            prev_PS2_power = PS2_power;
        end
    end                            
    
    % 显示最终结果
    disp('Results:');
    disp('Iter1 Iter2  PS2_fai(rad)  PS3_fai(rad)  PS5_fai(rad)  PS2_fai+up_fai0(rad)  PS3_fai+down_fai0(rad)  PS5_fai+IQ_fai0(rad)  PS2_Power  PS3_Power  PS5_Power');
    disp('------------------------------------------------------------------------------------------------------------------------------------------------------');
    for i = 1:size(results,1)
        row = results(i,:);
        fprintf('%4d %4d %14.4f %14.4f %14.4f %14.4f %22.4f %22.4f %10.6f %10.6f %10.6f\n', ...
                row(1), row(2), row(3), row(4), row(5), ...
                mod(row(3) + up_fai0,2*pi), mod(row(4) + down_fai0,2*pi), ...
                mod(row(5) + IQ_fai0,2*pi), row(6), row(7), row(8));
    end
end

% 辅助函数：相位扫描
function [best_phase, mid_power] = scan_phase_mid(PS2, PS3, PS5, ...
                                         IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                         I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 700); % 扫描1000个点
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
    
    % 可选：如果需要更精确的结果，可以在找到的点附近进行二次扫描
    if idx > 1 && idx < length(phase_points)
        fine_phase_points = linspace(phase_points(idx-1), phase_points(idx+1), 100);
        fine_power_values = zeros(size(fine_phase_points));
        
        for i = 1:length(fine_phase_points)
            phase = fine_phase_points(i);
            
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
            
            I = abs(layout_IQ(IQ_IN, R, IQ_fai0, ...
                 up_fai0, I_RF_phi, current_PS2, ...
                 down_fai0, Q_RF_phi, current_PS3, ...
                 current_PS5, k, b))^2;
            fine_power_values(i) = I;
        end
        
        [~, fine_idx] = min(abs(fine_power_values - target_power));
        best_phase = fine_phase_points(fine_idx);
        mid_power = fine_power_values(fine_idx);
    end
end

function [best_phase, min_power] = scan_phase_min(PS2, PS3, PS5, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 700); % 扫描100个点
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

function [best_phase, max_power] = scan_phase_max(PS2, PS3, PS5, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 700); % 扫描101个点
    max_power = 0;
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
        
        if current_power > max_power
            max_power = current_power;
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
