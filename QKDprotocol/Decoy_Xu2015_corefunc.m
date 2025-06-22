function [R_bitperpulse, R_bitpersecond, e_obs, e_x1_U, q, N, time] = Decoy_Xu2015_corefunc ... 
            (epsilon_cor, epsilon_sec, ...
             qZ, nZ, f, gate_width, width_3dB, ...
             k, pk, ...
             L, alpha, ita_Bob_detect, ...
             D_theta, e_mis_Z, f_EC, ...
             dc_count, deadtime)
    % 输入：
    % epsilon_cor, epsilon_sec: security criteria
    % qX: 选X基概率(efficient BB84), nX: 块长, f[Hz]: 频率, gate_width[s]: 门宽,
    % 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[m]: A到B距离, alpha[dB/m]: 光路损耗
    % ita_Bob_detect: Bob解码端总效率（假设通过PC使得探测器效率几乎一致）,包括滤波解码端超导等，不包括门宽
    % e_mis_Z: A端态制备的Z基错误率，用Z基对比度反算得到（文章里假设phi0是完美的), f_EC: 纠错效率
    % dc_count: 超导单通道每秒暗计数，是个1x3向量（三态协议，三个超导）, deadtime[s]: 死时间
    % 输出：
    % R_bitperpulse, R_bitpersecond:密钥率
    % e_obs: Z基比特误码率（用于成码的是Z基）
    % e_x1_U: 估计得到的用于成码的Z基的相位误码率
    % q: Z基和X基态的最大保真度
    % N: 为了获得nX长度的有限密钥，A需要发送的脉冲数
    % time: 为了获得nX长度的有限密钥，需要进行QKD的时间
    
    %% 定义光路到Bob解码的效率
    ita_ch = 10.^(-alpha*L/10);
    ita_sys = ita_ch*ita_Bob_detect; % 假设通过PC使得探测器效率几乎一致
     
    %% 基矢选择概率
    qX = 1 - qZ; % A选Z基概率
    Pz = qZ^2; % AB都选Z基概率
    Px = qX^2; % AB都选X基概率
    
    %% 暗计数概率，这玩意特别影响密钥率
    p_dc_Z = (dc_count(1) + dc_count(2))/2*gate_width/1;   % 假设超导前两个通道用来检偏Z  
    p_dc_X = dc_count(3)*gate_width/1;   % 假设超导第三个通道用来检偏X
    
    %% 卡门宽带来的效率
    [~, width_eff] = gaussian_pulse_analysis(width_3dB, gate_width);
    
    %% 计算B端收到每个强度态的概率（Z基）
    % 使用Finite-key analysis on the 1-decoy state QKD protocol的推导
    Dk_z = Pz*pk.*((1-exp(-ita_sys.*k))*width_eff + p_dc_Z);  % B端探测每个强度态的概率(不含后脉冲项)  
    Cdt_z =  1./(1+f*Dk_z*deadtime);    % B端探测器死时间修正
    Rk_z = Dk_z.*Cdt_z;  % B端探测每个强度态的概率(含后脉冲项)    

    %% 计算B端收到每个强度态的概率（X基）
    % 使用Finite-key analysis on the 1-decoy state QKD protocol的推导
    Dk_x = Px*pk.*((1-exp(-ita_sys.*k))*width_eff + p_dc_X);  % B端探测每个强度态的概率(不含后脉冲项)  
    Cdt_x =  1./(1+f*Dk_x*deadtime);    % B端探测器死时间修正
    Rk_x = Dk_x.*Cdt_x;  % B端探测每个强度态的概率(含后脉冲项)
    
    %% 计算A端制备的各个基上各个强度态数目
    N = nZ/sum(Rk_z);  % A端一共发送的态数目    
    n_Z = N.*Rk_z;   % B端接收的X基各强度态数目
      
    %% 利用制备态夹角计算mismatch基矢的态数目
    [~, delta_theta] = calculate_delta_theta(D_theta, [ita_Bob_detect, ita_Bob_detect], 10^(-10));
    [X_basis_state, Z_basis_state, p_0x0z, p_0x1z, p_0x0x, p_1x0z, p_1x1z, p_1x0x] ...
            = calculate_pdtec(delta_theta);
    q = compute_q(X_basis_state, Z_basis_state);    % 计算Z基和X基态的最大保真度
    
    n0x0zk = qZ*qX.*(p_dc_X/2 + p_0x0z*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    n0x1zk = qZ*qX.*(p_dc_X/2 + p_0x1z*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    n0x0xk = qX*qX.*(p_dc_X/2 + p_0x0x*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    
    n1x0zk = qZ*qX.*(p_dc_X/2 + p_1x0z*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    n1x1zk = qZ*qX.*(p_dc_X/2 + p_1x1z*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    n1x0xk = qX*qX.*(p_dc_X/2 + p_1x0x*(1-exp(-ita_ch *k))).*Cdt_x.*N;
    
    
    %% 计算Z基比特误码率
    e_k_z = Pz*pk.*(p_dc_Z/2 + e_mis_Z*(1-exp(-ita_ch *k))).*Cdt_z; % B端探测z基到且出错的概率
    
    m_Z = e_k_z.*N;  % B端接收的Z基各强度态出错的数目
    mZ = sum(m_Z); 
    
    e_obs = mZ/nZ;  
    
    %% 计算e_x1_U
    n0x0zk_plus = zeros(1,3);
    n0x0zk_minus = zeros(1,3);
    n0x1zk_plus = zeros(1,3);
    n0x1zk_minus = zeros(1,3);
    n0x0xk_plus = zeros(1,3);
    n0x0xk_minus = zeros(1,3);
    
    n1x0zk_plus = zeros(1,3);
    n1x0zk_minus = zeros(1,3);
    n1x1zk_plus = zeros(1,3);
    n1x1zk_minus = zeros(1,3);
    n1x0xk_plus = zeros(1,3);
    n1x0xk_minus = zeros(1,3);
    for i = 1:3
       n0x0zk_plus(i) = calculate_nzk_plus(k(i), pk(i), n0x0zk(i), sum(n0x0zk), epsilon_sec);
       n0x0zk_minus(i) = calculate_nzk_minus(k(i), pk(i), n0x0zk(i), sum(n0x0zk), epsilon_sec);
       n0x1zk_plus(i) = calculate_nzk_plus(k(i), pk(i), n0x1zk(i), sum(n0x1zk), epsilon_sec);
       n0x1zk_minus(i) = calculate_nzk_minus(k(i), pk(i), n0x1zk(i), sum(n0x1zk), epsilon_sec);
       n0x0xk_plus(i) = calculate_nxk_plus(k(i), pk(i), n0x0xk(i), sum(n0x0xk), epsilon_sec);
       n0x0xk_minus(i) = calculate_nxk_minus(k(i), pk(i), n0x0xk(i), sum(n0x0xk), epsilon_sec);

       n1x0zk_plus(i) = calculate_nzk_plus(k(i), pk(i), n1x0zk(i), sum(n1x0zk), epsilon_sec);
       n1x0zk_minus(i) = calculate_nzk_minus(k(i), pk(i), n1x0zk(i), sum(n1x0zk), epsilon_sec);
       n1x1zk_plus(i) = calculate_nzk_plus(k(i), pk(i), n1x1zk(i), sum(n1x1zk), epsilon_sec);
       n1x1zk_minus(i) = calculate_nzk_minus(k(i), pk(i), n1x1zk(i), sum(n1x1zk), epsilon_sec);
       n1x0xk_plus(i) = calculate_nxk_plus(k(i), pk(i), n1x0xk(i), sum(n1x0xk), epsilon_sec);
       n1x0xk_minus(i) = calculate_nxk_minus(k(i), pk(i), n1x0xk(i), sum(n1x0xk), epsilon_sec);
    end
    S_0x0z_1_plus = calculate_S_jxiz_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n0x0zk_plus(2), n0x0zk_minus(3));
    S_0x1z_1_plus = calculate_S_jxiz_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n0x1zk_plus(2), n0x1zk_minus(3));
    S_0x0x_1_plus = calculate_S_jx0x_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n0x0xk_plus(2), n0x0xk_minus(3));
    S_0x0z_1_minus =  calculate_S_jxiz_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n0x0zk_plus(1), n0x0zk_plus(2), n0x0zk_minus(2), n0x0zk_plus(3), n0x0zk_minus(3));
    S_0x1z_1_minus =  calculate_S_jxiz_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n0x1zk_plus(1), n0x1zk_plus(2), n0x1zk_minus(2), n0x1zk_plus(3), n0x1zk_minus(3));   
    S_0x0x_1_minus =  calculate_S_jx0x_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n0x0xk_plus(1), n0x0xk_plus(2), n0x0xk_minus(2), n0x0xk_plus(3), n0x0xk_minus(3));

    S_1x0z_1_plus = calculate_S_jxiz_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n1x0zk_plus(2), n1x0zk_minus(3));
    S_1x1z_1_plus = calculate_S_jxiz_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n1x1zk_plus(2), n1x1zk_minus(3));
    S_1x0x_1_plus = calculate_S_jx0x_1_plus(calculate_tau_n(k, pk, 1), k(2), k(3), n1x0xk_plus(2), n1x0xk_minus(3));
    S_1x0z_1_minus =  calculate_S_jxiz_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n1x0zk_plus(1), n1x0zk_plus(2), n1x0zk_minus(2), n1x0zk_plus(3), n1x0zk_minus(3));
    S_1x1z_1_minus =  calculate_S_jxiz_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n1x1zk_plus(1), n1x1zk_plus(2), n1x1zk_minus(2), n1x1zk_plus(3), n1x1zk_minus(3));   
    S_1x0x_1_minus =  calculate_S_jx0x_1_minus(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), ...
                            n1x0xk_plus(1), n1x0xk_plus(2), n1x0xk_minus(2), n1x0xk_plus(3), n1x0xk_minus(3));

    s_real_U = [S_0x0z_1_plus, S_0x1z_1_plus, S_0x0x_1_plus; ...
                S_1x0z_1_plus, S_1x1z_1_plus, S_1x0x_1_plus];
    s_real_L = [S_0x0z_1_minus, S_0x1z_1_minus, S_0x0x_1_minus; ...
                S_1x0z_1_minus, S_1x1z_1_minus, S_1x0x_1_minus];
            
    [s_vir_U, s_vir_L] = calculate_virtual_yields(qX, delta_theta, s_real_U, s_real_L);
    
    e_x1_U = calculate_phase_error_upper_bound(s_vir_U, s_vir_L);
    
    %% 计算SZ_0、SZ_1
    nzk_plus = zeros(1,3);
    nzk_minus = zeros(1,3);

    for i = 1:3
       nzk_plus(i) = calculate_nxk_plus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);
       nzk_minus(i) = calculate_nxk_minus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);       
    end
    SZ_0 = calculate_SZ_0(calculate_tau_n(k, pk, 0), k(2), k(3), nzk_plus(2), nzk_minus(3));
    SZ_1 = calculate_SZ_1(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), nzk_plus(1), nzk_minus(2), nzk_plus(3), SZ_0);
    
    %% 计算密钥率
    l = calculate_l(SZ_0, SZ_1, q, e_x1_U, f_EC*binary_entropy(e_obs)*nZ, epsilon_sec, epsilon_cor);
    R_bitperpulse= l/N;
    R_bitpersecond = l/N*f;
    time = N/f;
end

%% 
function [delta, delta_theta] = calculate_delta_theta(D_theta, ita, epsilon)
% 计算调制误差δ_θ的上界(式12)
%
% 输入参数:
%   D_theta - 一个nx2矩阵，其中D_theta(:,1)对应D_{1,θ}，D_theta(:,2)对应D_{2,θ}
%   ita - 探测器效率向量[η_{d1}, η_{d2}]
%   epsilon - 失败概率(如10^-10)
%   原文扣除暗计数没必要，假设计数足够大的话暗计数可以忽略不计

% 输出:
%   delta - 最大的delta_theta
%   delta_theta - 调制误差δ_θ的上界(弧度),是个1x2向量

    % 定义三态theta
    theta = [0; pi/2; pi];

    % 提取输入参数
    D1_theta = D_theta(:,1);  % D_{1,θ}的值
    D2_theta = D_theta(:,2);  % D_{2,θ}的值
    ita_d1 = ita(1);          % 探测器1的效率
    ita_d2 = ita(2);          % 探测器2的效率

    % 计算D_{1,0} (全局不对准)
    D1_0 = D1_theta(1);                        

    % 计算Δ项(使用Hoeffding不等式)
    Delta_D1_theta = sqrt(D1_theta/2 * log(1/epsilon));
    Delta_D2_theta = sqrt(D2_theta/2 * log(1/epsilon));
    Delta_D1_0 = sqrt(D1_0/2 * log(1/epsilon));


    % 计算分子和分母项
    numerator = ((D1_theta + Delta_D1_theta) - (D1_0 - Delta_D1_0)) / ita_d1;
    denominator = ((D2_theta - Delta_D2_theta) - (D1_0 + Delta_D1_0)) / ita_d2;


    % 计算arctan项
    atan_term = 2 * atan(sqrt(numerator ./ denominator));

    % 计算δ_θ (绝对值)
    delta_theta = abs(atan_term - theta);  % 这里假设预期相位θ=π
    delta_theta = delta_theta(2:3)';
    delta = max(delta_theta);
end

function [X_basis_state, Z_basis_state, p_0x0z, p_0x1z, p_0x0x, p_1x0z, p_1x1z, p_1x0x] ...
            = calculate_pdtec(delta_theta)
    % 提取 delta1 和 delta2
    delta1 = delta_theta(1);
    delta2 = delta_theta(2);
    
    % Pauli 矩阵
    sigma_z = [1, 0; 0, -1];
    I = eye(2);
    
    % 构造 rho0z = (I + sigma_z)/2
    rho0z = (I + sigma_z)/2;
    
    % 构造 rho1z = [sin^2(delta2), sin(delta2)cos(delta2); ...]
    rho1z = [sin(delta2)^2, sin(delta2)*cos(delta2);
             sin(delta2)*cos(delta2), cos(delta2)^2];
    
    % 构造 rho0x = 0.5*[1 + sin(2*delta1), cos(2*delta1); ...]
    rho0x = 0.5 * [1 + sin(2*delta1), cos(2*delta1);
                   cos(2*delta1), 1 - sin(2*delta1)];
    % 返回密度矩阵
    X_basis_state = {rho0x};
    Z_basis_state = {rho0z, rho1z};
    
    % 构造测量算子
    D_0x = 0.5 * [1, 1; 1, 1]; % |0_x⟩⟨0_x|
    D_1x = 0.5 * [1, -1; -1, 1]; % |1_x⟩⟨1_x|
    
    % 计算概率
    p_0x0z_raw = (1/6) * trace(D_0x * rho0z);
    p_0x1z_raw = (1/6) * trace(D_0x * rho1z);
    p_0x0x_raw = (1/6) * trace(D_0x * rho0x);
    
    p_1x0z_raw = (1/6) * trace(D_1x * rho0z);
    p_1x1z_raw = (1/6) * trace(D_1x * rho1z);
    p_1x0x_raw = (1/6) * trace(D_1x * rho0x);
    
    p_0x0z = p_0x0z_raw/(p_0x0z_raw + p_1x0z_raw);
    p_0x1z = p_0x1z_raw/(p_0x1z_raw + p_1x1z_raw);
    p_0x0x = p_0x0x_raw/(p_0x0x_raw + p_1x0x_raw);
    
    p_1x0z = p_1x0z_raw/(p_0x0z_raw + p_1x0z_raw);
    p_1x1z = p_1x1z_raw/(p_0x1z_raw + p_1x1z_raw);
    p_1x0x = p_1x0x_raw/(p_0x0x_raw + p_1x0x_raw);
end

function q = compute_q(X_basis, Z_basis)
    % 输入:
    % X_basis: X基矢的密度矩阵集合，cell数组，每个元素是一个2x2密度矩阵
    % Z_basis: Z基矢的密度矩阵集合，cell数组，每个元素是一个2x2密度矩阵
    %
    % 输出:
    % q: q = -log2(max |<psi_x|psi_z>|^2)

    max_fidelity = 0;

    % 遍历所有X基矢和Z基矢的组合
    for i = 1:length(X_basis)
        for j = 1:length(Z_basis)
            rho_x = X_basis{i};
            rho_z = Z_basis{j};

            % 计算保真度（假设是纯态）
            % 如果是纯态 |psi_x><psi_x| 和 |psi_z><psi_z|，保真度为 |<psi_x|psi_z>|^2
            % 提取态矢量（假设密度矩阵是纯态）
            [Vx, Dx] = eig(rho_x);
            [Vz, Dz] = eig(rho_z);
            psi_x = Vx(:, diag(Dx) == 1); % 提取本征值为1的态矢量
            psi_z = Vz(:, diag(Dz) == 1);

            % 计算内积的平方
            fidelity = abs(psi_x' * psi_z)^2;

            % 更新最大保真度
            if fidelity > max_fidelity
                max_fidelity = fidelity;
            end
        end
    end

    % 计算 q = -log2(max_fidelity)
    q = -log2(max_fidelity);
end

%%
function [sigma, energy_ratio] = gaussian_pulse_analysis(fwhm, gate_width)
% GAUSSIAN_PULSE_ANALYSIS 计算高斯脉冲的标准差和指定门宽内的能量占比
%
% 输入参数:
%   fwhm: 高斯脉冲的3dB全宽（FULL WIDTH AT HALF MAXIMUM）
%   gate_width: 以中心为基准的门的宽度（总宽度为 ±gate_width/2）
%
% 输出参数:
%   sigma: 高斯脉冲的标准差
%   energy_ratio: 在 gate_width 范围内的能量占比（0~1）

    % 1. 根据 FWHM 计算标准差 sigma
    % FWHM = 2 * sqrt(2 * log(2)) * sigma ≈ 2.355 * sigma
    sigma = fwhm / (2 * sqrt(2 * log(2)));  % 精确计算

    % 2. 计算门宽内的能量占比（使用误差函数 erf 计算积分）
    % 高斯函数在 [-a, a] 内的积分 = erf(a / (sqrt(2)*sigma))
    a = gate_width / 2;  % 半宽
    energy_ratio = erf(a / (sqrt(2) * sigma));
end

%%
function l = calculate_l(SZ_0, SZ_1, q, e_x1_U, lambda_EC, epsilon_sec, epsilon_cor)
    l = SZ_0 + SZ_1*q - SZ_1 * binary_entropy(e_x1_U) - lambda_EC - 6 * log2(26/epsilon_sec) - log2(2/epsilon_cor);
end

function e_x1_U = calculate_phase_error_upper_bound(s_vir_U, s_vir_L)
% 计算相位错误率上界 e_{x,1}^{U} (根据给定公式)
%
% 输入参数:
%   s_vir_U - 2x2矩阵，虚拟态概率上界:
%             s_vir_U(1,:) = [sz_0x0x_U, sz_0x1x_U] (s=0的情况)
%             s_vir_U(2,:) = [sz_1x0x_U, sz_1x1x_U] (s=1的情况)
%   s_vir_L - 2x2矩阵，虚拟态概率下界(格式同s_vir_U)
%
% 输出:
%   e_x1_U - 相位错误率上界 e_{x,1}^{U}

    % 从s_vir_U中提取分子所需项
    s_0x1x_vir_U = s_vir_U(1,2);  % s_{0_x|1_x,1}^{vir,U}
    s_1x0x_vir_U = s_vir_U(2,1);  % s_{1_x|0_x,1}^{vir,U}
    
    % 从s_vir_L中提取分母所需项
    s_0x0x_vir_L = s_vir_L(1,1);  % s_{0_x|0_x,1}^{vir,L}
    s_0x1x_vir_L = s_vir_L(1,2);  % s_{0_x|1_x,1}^{vir,L}
    s_1x0x_vir_L = s_vir_L(2,1);  % s_{1_x|0_x,1}^{vir,L}
    s_1x1x_vir_L = s_vir_L(2,2);  % s_{1_x|1_x,1}^{vir,L}
    
    % 计算分子和分母
    numerator = s_0x1x_vir_U + s_1x0x_vir_U;
    denominator = s_0x0x_vir_L + s_0x1x_vir_L + s_1x0x_vir_L + s_1x1x_vir_L;
    
    % 计算相位错误率上界
    e_x1_U = numerator / denominator;
end

function [s_vir_U, s_vir_L] = calculate_virtual_yields(qX, delta_theta, s_real_U, s_real_L)
% 计算虚拟态概率的上下界 (式8-9)
%
% 输入参数:
%   qX - Alice选择X基的概率 (qZ = 1-qX)
%   delta - 1x2向量 [δ1, δ2]，表示调制误差(弧度)
%   s_real_U - 2x3矩阵，包含实际上界测量值:
%              s_real_U(1,:) = [sz_0x0z_U, sz_0x1z_U, sx_0x0x_U] (s=0的情况)
%              s_real_U(2,:) = [sz_1x0z_U, sz_1x1z_U, sx_1x0x_U] (s=1的情况)
%   s_real_L - 2x3矩阵，包含实际下界测量值(格式同s_real_U)
%
% 输出:
%   s_vir_U - 2x2矩阵，虚拟态概率上界:
%             s_vir_U(1,:) = [sz_0x0x_U, sz_0x1x_U] (s=0的情况)
%             s_vir_U(2,:) = [sz_1x0x_U, sz_1x1x_U] (s=1的情况)
%   s_vir_L - 2x2矩阵，虚拟态概率下界(格式同s_vir_U)

    % 提取δ值和基选择概率
    delta1 = delta_theta(1);
    delta2 = delta_theta(2);
    qZ = 1 - qX;  % Alice选择Z基的概率

    % 构建矩阵A (式A8)
    A = [1, 1, 0;
         1, -cos(2*delta2), sin(2*delta2);
         1, sin(2*delta1), cos(2*delta1)];

    % 构建矩阵B (式A10)
    B = [(1+sin(delta2)), sin(delta2)*(1+sin(delta2)), cos(delta2)*(1+sin(delta2));
         (1-sin(delta2)), -sin(delta2)*(1-sin(delta2)), -cos(delta2)*(1-sin(delta2))] / 12;

    % 计算BA^{-1} (式A11)
    BA_inv = B / A;  % 等价于 B * inv(A)

    % 初始化输出矩阵
    s_vir_U = zeros(2,2);
    s_vir_L = zeros(2,2);

    % 对s=0和s=1分别计算上下界
    for s = 1:2
        % 上界计算 (式8)
        input_U = [2*qX*s_real_U(s,1); 
                  2*qX*s_real_U(s,2); 
                  qZ*s_real_U(s,3)];
        s_vir_s_U = BA_inv * input_U;
        s_vir_U(:,s) = s_vir_s_U ./ qZ;  % 除以qZ

        % 下界计算 (式9)
        input_L = [2*qX*s_real_L(s,1);
                  2*qX*s_real_U(s,2);  % 注意这里第二个元素用上界
                  qZ*s_real_U(s,3)];   % 第三个元素用上界
        s_vir_s_L = BA_inv * input_L;
        s_vir_L(:,s) = s_vir_s_L ./ qZ;  % 除以qZ
    end
end

%%
function S_jxiz_1_plus = calculate_S_jxiz_1_plus(tau_1, mu_2, mu_3, njxizk_mu2_plus, njxizk_mu3_minus)
    S_jxiz_1_plus = tau_1 * (njxizk_mu2_plus - njxizk_mu3_minus) / (mu_2 - mu_3);
end

function S_jxiz_1_minus = calculate_S_jxiz_1_minus(tau_0, tau_1, mu_1, mu_2, mu_3, ...
                            njxizk_mu1_plus, njxizk_mu2_plus, njxizk_mu2_minus, njxizk_mu3plus, njxizk_mu3_minus)
    S_jxiz_0 = calculate_SZ_0(tau_0, mu_2, mu_3, njxizk_mu2_plus, njxizk_mu3_minus);
    S_jxiz_1_minus = calculate_SZ_1(tau_0, tau_1, mu_1, mu_2, mu_3, njxizk_mu1_plus, njxizk_mu2_minus, njxizk_mu3plus, S_jxiz_0);
end

function S_jx0x_1_plus = calculate_S_jx0x_1_plus(tau_1, mu_2, mu_3, njxi0x_mu2_plus, njxi0x_mu3_minus)
    S_jx0x_1_plus = tau_1 * (njxi0x_mu2_plus - njxi0x_mu3_minus) / (mu_2 - mu_3);
end

function S_jx0x_1_minus = calculate_S_jx0x_1_minus(tau_0, tau_1, mu_1, mu_2, mu_3, ...
                            njxizk_mu1_plus, njx0xk_mu2_plus, njx0xk_mu2_minus, njx0xk_mu3plus, njx0xk_mu3_minus)
    S_jx0x_0 = calculate_SZ_0(tau_0, mu_2, mu_3, njx0xk_mu2_plus, njx0xk_mu3_minus);
    S_jx0x_1_minus = calculate_SZ_1(tau_0, tau_1, mu_1, mu_2, mu_3, njxizk_mu1_plus, njx0xk_mu2_minus, njx0xk_mu3plus, S_jx0x_0);
end

%%
function SZ_0 = calculate_SZ_0(tau_0, mu_2, mu_3, nZ_mu2_plus, nZ_mu3_minus)
    SZ_0 = tau_0 * (mu_2 * nZ_mu3_minus - mu_3 * nZ_mu2_plus) / (mu_2 - mu_3);
end

function SZ_1 = calculate_SZ_1(tau_0, tau_1, mu_1, mu_2, mu_3, n_plus_Z_mu1, n_minus_Z_mu2, n_plus_Z_mu3, SZ_0)
    SZ_1 = tau_1 * mu_1 * ...
          (n_minus_Z_mu2 - n_plus_Z_mu3 - ...
           ((mu_2^2 - mu_3^2) / mu_1^2) * (n_plus_Z_mu1 - SZ_0 / tau_0)) / ...
          (mu_1 * (mu_2 - mu_3) -  mu_2^2 + mu_3^2);
end

%%
function nxk_plus = calculate_nxk_plus(k, pk, nx_k, nx, epsilon_sec)
    nxk_plus = (exp(k)/pk) * (nx_k + sqrt(nx/2 * log(26/epsilon_sec)));
end

function nxk_minus = calculate_nxk_minus(k, pk, nx_k, nx, epsilon_sec)
    nxk_minus = (exp(k)/pk) * (nx_k - sqrt(nx/2 * log(26/epsilon_sec)));
end

function nzk_plus = calculate_nzk_plus(k, pk, nz_k, nz, epsilon_sec)
    nzk_plus = (exp(k)/pk) * (nz_k + sqrt(nz/2 * log(26/epsilon_sec)));
end

function nzk_minus = calculate_nzk_minus(k, pk, nz_k, nz, epsilon_sec)
    nzk_minus = (exp(k)/pk) * (nz_k - sqrt(nz/2 * log(26/epsilon_sec)));
end

%%
function tau_n = calculate_tau_n(K, Pk, n)
tau_n = 0;
    for i = 1:length(K)
        term = exp(-K(i)) * K(i)^n * Pk(i) / factorial(n);
        tau_n = tau_n + term;
    end
end

function entropy = binary_entropy(p)
    if p < 0 || p > 1
        entropy = 0;
    end
    entropy = -p * log2(p) - (1-p) * log2(1-p);
end
