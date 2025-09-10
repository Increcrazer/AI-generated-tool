function [R_bitperpulse, R_bitpersecond, e_obs, phi_X, nX] = Decoy_Lim2014_corefunc ... 
            (epsilon_cor, epsilon_sec, ...
             qX, N, f, gate_width, width_3dB, ...
             k, pk, ...
             L, alpha, eta_Bob_detect, ...
             e_mis_X, e_mis_Z, f_EC, ...
             p_ap, dc_count, deadtime)
    % 输入：
    % epsilon_cor, epsilon_sec: security criteria
    % qX: 选X基概率(efficient BB84), N: A发送的块长, f[Hz]: 频率, gate_width[s]: 门宽,
    % 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[km]: A到B距离, alpha[dB/km]: 光路损耗
    % eta_Bob_detect: Bob解码端总效率（假设通过PC使得探测器效率几乎一致）,包括滤波解码端超导等，不包括门宽
    % e_mis_X(Z): A端态制备的X(Z)基错误率，用X(Z)基对比度反算得到, f_EC: 纠错效率
    % p_ap: 后脉冲概率, dc_count[/s]: 超导单通道暗计数, 是个1x4向量（四个超导）,deadtime[s]: 死时间
    % 输出：
    % R_bitperpulse, R_bitpersecond:密钥率
    % e_obs: X基比特误码率（用于成码的是X基）
    % phi_X: 估计得到的用于成码的X基的相位误码率
    % nX: B接收到的X基块长
    
    %% 定义光路到Bob解码的效率
    eta_ch = 10.^(-alpha.*L./10);
    [~, width_eff] = gaussian_pulse_analysis(width_3dB, gate_width); % 卡门宽带来的效率
    eta_sys = eta_ch.*eta_Bob_detect.*width_eff; % 假设通过PC使得探测器效率几乎一致

    %% 暗计数概率，这玩意特别影响密钥率 
    p_dc = calculate_p_dc(dc_count, gate_width); 
    p_dc_X = p_dc(1);
    p_dc_Z = p_dc(2);
   
    %% 基矢选择概率
    qZ = 1 - qX; % 选Z基概率
    Px = qX.^2; % AB都选X基概率
    Pz = qZ.^2; % AB都选Z基概率
    
    %% 计算B端收到每个强度态的概率（X基）
%     %%% Lim 2014 %%%
%     Dk_x = Px.*pk.*(1 - (1 - 2.*p_dc_X).*exp(-eta_sys.*k));
%     Rk_x = Dk_x.*(1 + p_ap);
    
    %%% Rusca 2018 %%%
    Dk_x = Px.*pk.*(1 - (1 - 2.*p_dc_X).*exp(-eta_sys.*k));
    Cdt_x =  1./(1 + f.*sum(Dk_x).*deadtime);    % B端探测器死时间修正
    Rk_x = Cdt_x.*Dk_x;  % B端探测每个强度态的概率(含死时间修正)
    
    %% 计算B端收到每个强度态的概率（Z基）
%     %%% Lim 2014 %%%
%     Dk_z = Pz.*pk.*(1 - (1 - 2.*p_dc_Z).*exp(-eta_sys.*k));
%     Rk_z = Dk_z.*(1 + p_ap);
    
    %%% Rusca 2018 %%%
    Dk_z = Pz.*pk.*(1 - (1 - 2.*p_dc_Z).*exp(-eta_sys.*k));
    Cdt_z =  1./(1 + f.*sum(Dk_z).*deadtime);    % B端探测器死时间修正
    Rk_z = Cdt_z.*Dk_z;  % B端探测每个强度态的概率(含死时间修正)
    
    %% 计算计算B端各强度态数目
    n_X = N.*Rk_x;   % B端接收的X基各强度态数目
    n_Z = N.*Rk_z;   % B端接收的Z基各强度态数目
    nX = sum(n_X);  % B端一共收到的X基态数目
    nZ = sum(n_Z);  % B端一共收到的Z基态数目
    
    %% 计算B端错误率
%     %%% Lim 2014 %%%
%     ek_x = Px.*pk.*(p_dc_X + e_mis_X.*(1-exp(-eta_sys .*k)) + p_ap.*Dk_x./2); % B端探测x基到且出错的概率   
%     ek_z = Pz.*pk.*(p_dc_Z + e_mis_Z.*(1-exp(-eta_sys .*k)) + p_ap.*Dk_z./2); % B端探测z基到且出错的概率
    
    %%% Rusca 2018 %%%
    ek_x = Px.*pk.*(p_dc_X + e_mis_X.*(1-exp(-eta_sys .*k))).*Cdt_x; % B端探测x基到且出错的概率   
    ek_z = Pz.*pk.*(p_dc_Z + e_mis_Z.*(1-exp(-eta_sys .*k))).*Cdt_z; % B端探测z基到且出错的概率
    
    %% 计算计算B端各强度态误码数目
    m_X = N.*ek_x;  % B端接收的X基各强度态出错的数目
    m_Z = N.*ek_z;  % B端接收的Z基各强度态出错的数目
    mX = sum(m_X);
    mZ = sum(m_Z); 

    e_obs = mX./nX;

    %% 计算phi_X、SX_0、SX_1    
    nxk_plus = zeros(1,3);
    nxk_minus = zeros(1,3);
    nzk_plus = zeros(1,3);
    nzk_minus = zeros(1,3);
    mzk_plus = zeros(1,3);
    mzk_minus = zeros(1,3);

    for i = 1:3
       nxk_plus(i) =  calculate_nxk_plus(k(i), pk(i), n_X(i), nX, epsilon_sec);
       nxk_minus(i) = calculate_nxk_minus(k(i), pk(i), n_X(i), nX, epsilon_sec);
       nzk_plus(i) = calculate_nxk_plus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);
       nzk_minus(i) = calculate_nxk_minus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);       
       mzk_plus(i) = calculate_mzk_plus(k(i), pk(i), m_Z(i), mZ, epsilon_sec);
       mzk_minus(i) = calculate_mzk_minus(k(i), pk(i), m_Z(i), mZ, epsilon_sec);
    end
    SX_0 = calculate_SX_0(calculate_tau_n(k, pk, 0), k(2), k(3), nxk_plus(2), nxk_minus(3));
    SX_1 = calculate_SX_1(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), nxk_plus(1), nxk_minus(2), nxk_plus(3), SX_0);

    SZ_0 = calculate_SX_0(calculate_tau_n(k, pk, 0), k(2), k(3), nzk_plus(2), nzk_minus(3));
    SZ_1 = calculate_SX_1(calculate_tau_n(k, pk, 0), calculate_tau_n(k, pk, 1), k(1), k(2), k(3), nzk_plus(1), nzk_minus(2), nzk_plus(3), SZ_0);
             
    nuZ_1 = calculate_nu_Z_1(calculate_tau_n(k, pk, 1), mzk_plus(2), mzk_minus(3), k(2), k(3));
    phi_X = calculate_phi_X(SX_1, nuZ_1, SZ_1, epsilon_sec);
    
    %% 计算密钥率
    if isfinite(phi_X) && phi_X >= 0 && phi_X <= 0.5 && e_obs >= 0 && e_obs <= 0.5 % 添加有效性检查
        l = real(calculate_l(SX_0, SX_1, phi_X, f_EC.*binary_entropy(e_obs).*nX, epsilon_sec, epsilon_cor));
    else
        l = 0;
    end

    if l < 0 || l >= N
        l = 0;
    end
    
    R_bitperpulse= l./N;
    R_bitpersecond = l./N.*f;
end

%%
function p_dc = calculate_p_dc(dc_count, gate_width)
    p_dc_X = (dc_count(1) + dc_count(2))./2.*gate_width./1;   % 假设超导前两个通道用来检偏X  
    p_dc_Z = (dc_count(3) + dc_count(4))./2.*gate_width./1;   % 假设超导后两个通道用来检偏Z 
    p_dc = [p_dc_X, p_dc_Z];
end

%%
function [sigma, energy_ratio] = gaussian_pulse_analysis(fwhm, gate_width)
% GAUSSIAN_PULSE_ANALYSIS 计算高斯脉冲的标准差和指定门宽内的能量占比
%
% 输入参数:
%   fwhm: 高斯脉冲的3dB全宽（FULL WIDTH AT HALF MAXIMUM）
%   gate_width: 以中心为基准的门的宽度（总宽度为 ±gate_width./2）
%
% 输出参数:
%   sigma: 高斯脉冲的标准差
%   energy_ratio: 在 gate_width 范围内的能量占比（0~1）

    % 1. 根据 FWHM 计算标准差 sigma
    % FWHM = 2 .* sqrt(2 .* log(2)) .* sigma ≈ 2.355 .* sigma
    sigma = fwhm ./ (2 .* sqrt(2 .* log(2)));  % 精确计算

    % 2. 计算门宽内的能量占比（使用误差函数 erf 计算积分）
    % 高斯函数在 [-a, a] 内的积分 = erf(a ./ (sqrt(2).*sigma))
    a = gate_width ./ 2;  % 半宽
    energy_ratio = erf(a ./ (sqrt(2) .* sigma));
end

%%
function l = calculate_l(SX_0, SX_1, phiX, lambda_EC, epsilon_sec, epsilon_cor)
    l = SX_0 + SX_1 - SX_1 .* binary_entropy(phiX) - lambda_EC - 6 .* log2(21./epsilon_sec) - log2(2./epsilon_cor);
end

function SX_0 = calculate_SX_0(tau_0, mu_2, mu_3, nX_mu2_plus, nX_mu3_minus)
    SX_0 = tau_0 .* (mu_2 .* nX_mu3_minus - mu_3 .* nX_mu2_plus) ./ (mu_2 - mu_3);
end

function SX_1 = calculate_SX_1(tau_0, tau_1, mu_1, mu_2, mu_3, n_plus_X_mu1, n_minus_X_mu2, n_plus_X_mu3, SX_0)
    SX_1 = tau_1 .* mu_1 .* ...
          (n_minus_X_mu2 - n_plus_X_mu3 - ...
           ((mu_2.^2 - mu_3.^2) ./ mu_1.^2) .* (n_plus_X_mu1 - SX_0 ./ tau_0)) ./ ...
          (mu_1 .* (mu_2 - mu_3) -  mu_2.^2 + mu_3.^2);
end

function phi_X = calculate_phi_X(SX_1, nuZ_1, SZ_1, epsilon_sec)
    phi_X = (nuZ_1 ./ SZ_1) + calculate_gamma(epsilon_sec, nuZ_1./SZ_1, SZ_1, SX_1);
end

function nu_Z_1 = calculate_nu_Z_1(tau1, m_Z_mu2, m_Z_mu3, mu2, mu3)
    nu_Z_1 = tau1 .* (m_Z_mu2 - m_Z_mu3) ./ (mu2 - mu3);
end

%%
function nxk_plus = calculate_nxk_plus(k, pk, nx_k, nx, epsilon_sec)
    nxk_plus = (exp(k)./pk) .* (nx_k + sqrt(nx./2 .* log(21./epsilon_sec)));
end

function nxk_minus = calculate_nxk_minus(k, pk, nx_k, nx, epsilon_sec)
    nxk_minus = (exp(k)./pk) .* (nx_k - sqrt(nx./2 .* log(21./epsilon_sec)));
end

%%
function mzk_plus = calculate_mzk_plus(k, pk, mz_k, mz, epsilon_sec)
    mzk_plus = (exp(k)./pk) .* (mz_k + sqrt((mz./2) .* log(21./epsilon_sec)));
end

function mzk_minus = calculate_mzk_minus(k, pk, mz_k, mz, epsilon_sec)
    mzk_minus = (exp(k)./pk) .* (mz_k - sqrt((mz./2) .* log(21./epsilon_sec)));
end

%%
function gamma = calculate_gamma(a, b, c, d)
    eqa = ((c + d) .* (1 - b) .* b) ./ (c .* d .* log10(2)) .* log2(((c + d) .* 21.^2) ./ ((c .* d .* (1 - b) .* b) .* a.^2));
    if eqa >= 0
        gamma = sqrt(eqa);
    else
        gamma = inf;
    end
end

function tau_n = calculate_tau_n(K, Pk, n)
tau_n = 0;
    for i = 1:length(K)
        term = exp(-K(i)) .* K(i).^n .* Pk(i) ./ factorial(n);
        tau_n = tau_n + term;
    end
end

function entropy = binary_entropy(p)
    if p < 0 || p > 1
        entropy = 0;
    end
    entropy = -p .* log2(p) - (1-p) .* log2(1-p);
end
