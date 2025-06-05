function find_IQ_phase_method()
    % 一种找到I和Q消光点的迭代算法
    % step 1收敛有可能出问题
    % 加上step 2精度更高
    
    % step 1可能只是让I路和Q路相位相同，有待研究
    
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
    num_iterations = 20;

    % 存储结果
    results = zeros(num_iterations, 4); % PS2, PS3, PS5, Power

    % step1:优化循环
    for iter = 1:num_iterations       
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
       
        % 扫描PS5_fai
        [PS5_fai, PS5_max_power] = scan_phase_max(PS2_fai, PS3_fai, PS5_fai, ...
                                         IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                         I_RF_phi, Q_RF_phi, k, b, 'PS5');
        PS5_fai = mod(PS5_fai,2*pi);        
        % 存储结果
        results(iter,:) = [PS2_fai, PS3_fai, PS5_fai, PS5_max_power];
    end
    % PS3设置为step1消光，寻找PS2通光点
    [PS2_fai_open, ~] = scan_phase_max(PS2_fai, PS3_fai, PS5_fai, ...
                                     IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                     I_RF_phi, Q_RF_phi, k, b, 'PS2');
    % PS2设置为step1消光，寻找PS3通光点
    [PS3_fai_open, ~] = scan_phase_max(PS2_fai, PS3_fai, PS5_fai, ...
                                     IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                     I_RF_phi, Q_RF_phi, k, b, 'PS3');  
    % PS2、PS3设置为通光点，寻找PS5通光点
    [PS5_fai_open, ~] = scan_phase_max(PS2_fai_open, PS3_fai_open, PS5_fai, ...
                                     IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                     I_RF_phi, Q_RF_phi, k, b, 'PS5');  
    % PS5、PS3设置为通光点，寻找PS2消光点
    [PS2_fai_close, ~] = scan_phase_min(PS2_fai, PS3_fai_open, PS5_fai_open, ...
                                     IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                     I_RF_phi, Q_RF_phi, k, b, 'PS2');    
    % PS5、PS2设置为通光点，寻找PS3消光点
    [PS3_fai_close, ~] = scan_phase_min(PS2_fai_open, PS3_fai, PS5_fai_open, ...
                                 IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                 I_RF_phi, Q_RF_phi, k, b, 'PS3');    
    
    % 显示最终结果
    disp('Step1 Results:');
    disp('Iteration  PS2_fai(rad)  PS3_fai(rad)  PS5_fai(rad)  PS2_fai+up_fai0(rad)  PS3_fai+down_fai0(rad)  PS5_fai+IQ_fai0(rad)  Power');
    disp('------------------------------------------------------------------------------------------------------------------------------');
    for i = 1:size(results,1)
        fprintf('%4d %14.4f %14.4f %14.4f %14.4f %22.4f %22.4f %16.6f\n', ...
                i, results(i,1), results(i,2), results(i,3), mod(results(i,1) + up_fai0,2*pi), ...
            mod(results(i,2) + down_fai0,2*pi), mod(results(i,3) + IQ_fai0,2*pi), results(i,4));
    end
    disp('------------------------------------------------------------------------------------------------------------------------------');
    disp('Step2 Results:');
    disp('------------------------------------------------------------------------');
    disp('PS2_fai(rad)  PS3_fai(rad)  PS2_fai+up_fai0(rad)  PS3_fai+down_fai0(rad)');
    fprintf('%8.4f %14.4f %20.4f %20.4f\n', ...
        mod(PS2_fai_close,2*pi), mod(PS3_fai_close,2*pi), mod(PS2_fai_close + up_fai0,2*pi), mod(PS3_fai_close + down_fai0,2*pi));
    disp('------------------------------------------------------------------------');
end

% 辅助函数：相位扫描
function [best_phase, min_power] = scan_phase_min(PS2, PS3, PS5, ...
                                             IQ_IN, R, IQ_fai0, up_fai0, down_fai0, ...
                                             I_RF_phi, Q_RF_phi, k, b, target)
    phase_points = linspace(0, 2*pi, 1000); % 扫描100个点
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
    phase_points = linspace(0, 2*pi, 101); % 扫描101个点
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
    0,PS5_fai, ...
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
