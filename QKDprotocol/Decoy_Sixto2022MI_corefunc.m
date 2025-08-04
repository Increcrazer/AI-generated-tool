function R_bitperpulse = Decoy_Sixto2022MI_corefunc ... 
            (qZ, ncut, gate_width, width_3dB, ...
             k, pk, ...
             L, alpha, eta_Bob_detect, ...
             delta_A, f_EC, ...
             pd, delta_max)
    % 输入：
    % qZ: 选Z基概率(efficient BB84), ncut: 诱骗态截止光子数, gate_width[s]: 门宽, 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[m]: A到B距离, alpha[dB/m]: 光路损耗
    % eta_Bob_detect: Bob解码端总效率（假设通过PC使得探测器效率几乎一致）,包括滤波解码端超导等，不包括门宽
    % delta_A: 偏振基偏差角度, f_EC: 纠错效率
    % pd: 暗计数概率, delta_max: 最大偏差
    % 输出：
    % R_bitperpulse:密钥率
    
    %%
    qX = 1- qZ; 
    n = 0:ncut;
    
    %% 计算每个真实强度的bound （假设相对便宜参数相等，由delta_max决定）
    k_plus = k.*(1 + delta_max);
    k_minus = k.*(1 - delta_max);
         
    %% 定义光路到Bob解码的效率
    eta_ch = 10.^(-alpha.*L./10);
    [~, width_eff] = gaussian_pulse_analysis(width_3dB, gate_width); % 卡门宽带来的效率
    eta_sys = eta_ch.*eta_Bob_detect.*width_eff; % 假设通过PC使得探测器效率几乎一致

    %% 参数计算部分
    [Dk_B, ek_B] = calculate_De_B(eta_sys, pd, k, delta_A);
    [y_tilde, h_tilde] = calculate_hy_tilde(delta_A, eta_sys, pd, n);

    CStao = caculate_CStao(pk, k_plus, k_minus);
    [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
            = caculate_CScmts(CStao, y_tilde, h_tilde, n);
    
    %% 线性规划部分
    [~, ZS1_bar_low] = linprog_sigma_y1h1(qZ, pk, k_plus, k_minus, Dk_B, n, c_plus, c_minus, m_plus, m_minus);
    [~, XS1_bar_low] = linprog_sigma_y1h1(qX, pk, k_plus, k_minus, Dk_B, n, c_plus, c_minus, m_plus, m_minus);
    [~, EXS1_bar_up] = linprog_sigma_h1h1(qX, pk, k_plus, k_minus, ek_B, n, t_plus, t_minus, s_plus, s_minus);

    %% 计算汇总值
    [Zkk_bar, ~] = caculate_ZXkk_bar(qZ, qX, pk, Dk_B);
    ZS_bar = sum(Zkk_bar(1,:));
    [EZkk_bar, ~] = caculate_Ekk_bar(qZ, qX, pk, ek_B);
    EZS_bar = sum(EZkk_bar(1,:));

    %% 计算密钥率
    R_bitperpulse  = caculate_SKR(ZS1_bar_low, XS1_bar_low, EXS1_bar_up, ZS_bar, EZS_bar./ZS_bar, f_EC);
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
function infinite_key = caculate_SKR(ZS1_bar_low, XS1_bar_low, EXS1_bar_up, ZS_bar, Etol, f_EC)
    infinite_key = ZS1_bar_low.*(1 - binary_entropy(EXS1_bar_up./XS1_bar_low)) - f_EC.*ZS_bar.*binary_entropy(Etol);
end


%%
function [Zkk_bar, Xkk_bar] = caculate_ZXkk_bar(qZ, qX, pk, Dk_B)
    Zkk_bar = zeros(3,3);
    Xkk_bar = zeros(3,3);
    for a = 1:3
        for c = 1:3
            Zkk_bar(a,c) = qZ.^2.*pk(a).*pk(c).*Dk_B(a);
            Xkk_bar(a,c) = qX.^2.*pk(a).*pk(c).*Dk_B(a);
        end
    end   
end

function [EZkk_bar, EXkk_bar] = caculate_Ekk_bar(qZ, qX, pk, ek_B)
    EZkk_bar = zeros(3,3);
    EXkk_bar = zeros(3,3);
    for a = 1:3
        for c = 1:3
            EZkk_bar(a,c) = qZ.^2.*pk(a).*pk(c).*ek_B(a);
            EXkk_bar(a,c) = qX.^2.*pk(a).*pk(c).*ek_B(a);
        end
    end   
end

%%
function [nu_opt, fval] = linprog_sigma_y1h1(qZ, pk, k_plus, k_minus, Dk_B, n, c_plus, c_minus, m_plus, m_minus)
    % 高精度QKD线性规划问题求解
    
    % 目标函数构建
    f = zeros(3*3*length(n), 1);
    for i = 1:3
        f((i-1)*length(n)+2) = qZ^2 * pk(1) * pk(i) * k_minus(1) * exp(-k_minus(1));
    end
    
    % 使用更高精度的约束构建方式
    total_vars = 3*3*length(n);
    
    % 预计算约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    
    % 预分配约束矩阵
    A = zeros(num_constr1 + num_constr2 + num_constr3 + num_constr4, total_vars);
    b = zeros(num_constr1 + num_constr2 + num_constr3 + num_constr4, 1);
    
    % 约束1: Dk_B约束
    constr_idx = 1;
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pos = (a-1)*3*length(n) + (c-1)*length(n);
            row(pos+1) = exp(-k_plus(a));
            for i = 2:length(n)
                ni = n(i);
                fact_ni = factorial(ni);
                row(pos+i) = exp(-k_minus(a)) * k_minus(a)^ni / fact_ni;
            end
            A(constr_idx,:) = row;
            b(constr_idx) = Dk_B(a);
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束2: 1-Dk_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pos = (a-1)*3*length(n) + (c-1)*length(n);
            row(pos+1) = -exp(-k_minus(a));

            % 预计算求和部分
            sum_terms = 0;
            for i = 2:length(n)
                ni = n(i);
                % 使用对数变换计算每一项
                log_term = -k_plus(a) + ni*log(k_plus(a)) - sum(log(1:ni));
                term = exp(log_term);
                sum_terms = sum_terms + term;
                row(pos+i) = -term;  % 直接使用计算好的term
            end

            A(constr_idx,:) = row;
            b(constr_idx) = 1 - exp(-k_plus(a)) - Dk_B(a) - sum_terms;
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束3: m_plus约束
    for a = 1:3
        for c = 1:3
            for bb = 1:3
                if bb == a, continue; end
                for i = 1:length(n)
                    row = zeros(1, total_vars);
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n);
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n);
                    row(pos_a+i) = -m_plus(a,bb,c,i);
                    row(pos_bb+i) = 1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = c_plus(a,bb,c,i);
                    constr_idx = constr_idx + 1;
                end
            end
        end
    end
    
    % 约束4: m_minus约束
    for a = 1:3
        for c = 1:3
            for bb = 1:3
                if bb == a, continue; end
                for i = 1:length(n)
                    row = zeros(1, total_vars);
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n);
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n);
                    row(pos_a+i) = m_minus(a,bb,c,i);
                    row(pos_bb+i) = -1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = -c_minus(a,bb,c,i);
                    constr_idx = constr_idx + 1;
                end
            end
        end
    end
    
    % 高精度求解选项
    options = optimoptions('linprog', ...
        'Algorithm', 'dual-simplex', ...  % 保持对偶单纯形法
        'ConstraintTolerance', 1e-9, ...  % 约束容差从默认1e-4提高到1e-8
        'OptimalityTolerance', 1e-9, ...  % 最优性容差从默认1e-6提高到1e-8
        'MaxIterations', 100000, ...      % 最大迭代次数提高到10万次
        'Display', 'final');               % 显示迭代过程便于调试
    
    % 求解线性规划
    [nu_opt, fval] = linprog(f, A, b, [], [], zeros(total_vars,1), ones(total_vars,1), [], options);
end

%%
function [nu_opt, fval] = linprog_sigma_h1h1(qX, pk, k_plus, k_minus, ek_B, n, t_plus, t_minus, s_plus, s_minus)
    % 高精度QKD线性规划问题求解
    
    % 目标函数构建
    f = zeros(3*3*length(n), 1);
    for i = 1:3
        f((i-1)*length(n)+2) = qX^2 * pk(1) * pk(i) * k_plus(1) * exp(-k_plus(1));
    end
    
    % 使用更高精度的约束构建方式
    total_vars = 3*3*length(n);
    
    % 预计算约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    
    % 预分配约束矩阵
    A = zeros(num_constr1 + num_constr2 + num_constr3 + num_constr4, total_vars);
    b = zeros(num_constr1 + num_constr2 + num_constr3 + num_constr4, 1);
    
    % 约束1: ek_B约束
    constr_idx = 1;
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pos = (a-1)*3*length(n) + (c-1)*length(n);
            row(pos+1) = exp(-k_plus(a));
            for i = 2:length(n)
                ni = n(i);
                fact_ni = factorial(ni);
                row(pos+i) = exp(-k_minus(a)) * k_minus(a)^ni / fact_ni;
            end
            A(constr_idx,:) = row;
            b(constr_idx) = ek_B(a);
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束2: 1-ek_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pos = (a-1)*3*length(n) + (c-1)*length(n);
            row(pos+1) = -exp(-k_minus(a));

            % 计算sum(exp(-k_plus(a)).*k_plus(a).^n(2:end)./factorial(n(2:end)))
            sum_terms = 0;
            for i = 2:length(n)
                ni = n(i);
                % 使用对数变换计算每一项，避免大数阶乘
                log_term = -k_plus(a) + ni*log(k_plus(a)) - sum(log(1:ni));
                term = exp(log_term);
                sum_terms = sum_terms + term;
                row(pos+i) = -term;  % 直接使用计算好的term
            end

            A(constr_idx,:) = row;
            b(constr_idx) = 1 - exp(-k_plus(a)) - ek_B(a) - sum_terms;
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束3: s_plus约束
    for a = 1:3
        for c = 1:3
            for bb = 1:3
                if bb == a, continue; end
                for i = 1:length(n)
                    row = zeros(1, total_vars);
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n);
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n);
                    row(pos_a+i) = -s_plus(a,bb,c,i);
                    row(pos_bb+i) = 1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = t_plus(a,bb,c,i);
                    constr_idx = constr_idx + 1;
                end
            end
        end
    end
    
    % 约束4: s_minus约束
    for a = 1:3
        for c = 1:3
            for bb = 1:3
                if bb == a, continue; end
                for i = 1:length(n)
                    row = zeros(1, total_vars);
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n);
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n);
                    row(pos_a+i) = s_minus(a,bb,c,i);
                    row(pos_bb+i) = -1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = -t_minus(a,bb,c,i);
                    constr_idx = constr_idx + 1;
                end
            end
        end
    end
    
    % 高精度求解选项
    options = optimoptions('linprog', ...
        'Algorithm', 'dual-simplex', ...  % 保持对偶单纯形法
        'ConstraintTolerance', 1e-9, ...  % 约束容差从默认1e-4提高到1e-8
        'OptimalityTolerance', 1e-9, ...  % 最优性容差从默认1e-6提高到1e-8
        'MaxIterations', 100000, ...      % 最大迭代次数提高到10万次
        'Display', 'final');               % 显示迭代过程便于调试
    
    % 求解线性规划
    [nu_opt, fval] = linprog(-f, A, b, [], [], zeros(total_vars,1), ones(total_vars,1), [], options);
    fval = -fval;
end

%%
function [Dk_B, ek_B] = calculate_De_B(eta, pd, k, delta_A)
    % 计算探测率Dk_B和错误率ek_B（支持k为向量输入）
    % 输入参数：
    %   eta     : 系统总效率 (η = η_def .* η_ch)
    %   pd      : 探测器暗计数概率
    %   k      : 强度向量
    %   delta_A : 偏振失配角 (rad)
    % 输出：
    %   Dk_B    : B端探测率向量 1x3
    %   ek_B    : B端错误率向量 1x3

    % 计算探测率 Dk_B (公式类似文中y_tilde)
    Dk_B = 1 - (1 - pd).^2 .* exp(-eta .* k);
    
    % 计算中间参数 h (公式37)
    h = (exp(-eta .* k .* cos(delta_A).^2) - exp(-eta .* k .* sin(delta_A).^2)) ./ 2;
    
    % 计算错误率 ek_B (根据公式36推导)
    term1 = pd.^2 ./ 2;  % 暗计数同时发生的项
    term2 = pd .* (1 - pd) .* (1 + h);  % 单边暗计数项
    term3 = (1 - pd).^2 .* (0.5 + h - 0.5 .* exp(-eta .* k));  % 无暗计数项
    
    ek_B = term1 + term2 + term3;
end

function [y_tilde, h_tilde] = calculate_hy_tilde(delta_A, eta, pd, n)
    % 计算量子密钥分发中的产出率和错误率参考值
    % 输入参数：
    %   delta_A : 偏振失配角 (rad)
    %   eta     : 系统总效率 (包括信道损耗和探测器效率)
    %   pd      : 探测器暗计数概率
    %   n       : 光子数
    % 输出：
    %   y_tilde : 产出率 \tilde{y}_{n,a,c} (与a,c无关) 单标量
    %   h_tilde : 错误率 \tilde{h}_{n,a,c} (与a,c无关) 单标量

    % 计算基础概率 (公式CI)
    p00 = (1 - eta).^n;
    p01 = (eta .* sin(delta_A).^2 + 1 - eta).^n - p00;
    p10 = (eta .* cos(delta_A).^2 + 1 - eta).^n - p00;
    p11 = 1 - p00 - p01 - p10;

    % 计算不同暗计数条件下的错误概率 (公式C2)
    p_err_A = p01 + 0.5 .* p11;
    p_err_B = 0.5 .* (p01 + p11);
    p_err_C = p00 + p01 + 0.5 .* (p10 + p11);
    p_err_D = 0.5;

    % 计算错误率 (公式C3)
    h_tilde = (1 - pd).^2 .* p_err_A + ...
              pd .* (1 - pd) .* (p_err_B + p_err_C) + ...
              pd.^2 .* p_err_D;

    % 计算产出率 (公式C4)
    y_tilde = 1 - (1 - pd).^2 .* p00;
end

function CStao = caculate_CStao(pk, k_plus, k_minus) 
    % CStao 3x3x3
    CStao = zeros(3,3,3);
        for i = 1:3 % tao_a
            for j = 1:3 % tao_b
                for k =1:3  % tao_c
                    sigma = 0;
                    for r = 1:3
                        sigma = sigma + pk(r).*(exp(-k_minus(r)) - exp(-k_plus(r)));
                    end
                    CStao(i,j,k) = (1 - sigma).^2;
                end
            end
        end
end

function [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
        = caculate_CScmts(CStao, y_tilde, h_tilde, n)
    % cmts 3x3x3x(ncut+1)
    c_plus = zeros(3,3,3,length(n));
    c_minus = zeros(3,3,3,length(n));
    m_plus = zeros(3,3,3,length(n));
    m_minus = zeros(3,3,3,length(n));
    t_plus = zeros(3,3,3,length(n));
    t_minus = zeros(3,3,3,length(n));
    s_plus = zeros(3,3,3,length(n));
    s_minus = zeros(3,3,3,length(n));
    for i = 1:3 % cmts_a
        for j = 1:3 % cmts_b
            for k = 1:3  % cmts_c
                for l = 1:length(n)
                    c_plus(i,j,k,l) = G_plus(y_tilde(l),CStao(i,j,k)) - Gprime_plus(y_tilde(l),CStao(i,j,k)).*y_tilde(l);
                    c_minus(i,j,k,l) = G_minus(y_tilde(l),CStao(i,j,k)) - Gprime_minus(y_tilde(l),CStao(i,j,k)).*y_tilde(l);
                    m_plus(i,j,k,l) = Gprime_plus(y_tilde(l),CStao(i,j,k));
                    m_minus(i,j,k,l) = Gprime_minus(y_tilde(l),CStao(i,j,k));
                    t_plus(i,j,k,l) = G_plus(h_tilde(l),CStao(i,j,k)) - Gprime_plus(h_tilde(l),CStao(i,j,k)).*h_tilde(l);
                    t_minus(i,j,k,l) = G_minus(h_tilde(l),CStao(i,j,k)) - Gprime_minus(h_tilde(l),CStao(i,j,k)).*h_tilde(l);
                    s_plus(i,j,k,l) = Gprime_plus(h_tilde(l),CStao(i,j,k));
                    s_minus(i,j,k,l) = Gprime_minus(h_tilde(l),CStao(i,j,k));                    
                end
            end
        end
    end
end

%%
function G = Gprime_minus(y, z)
    % 计算 G'_-(y,z)，输入 y 和 z 必须为标量
    if y > 1 - z
        term = sqrt(z .* (1 - z) ./ (y .* (1 - y)));
        G = -1 + 2.*z - (1 - 2.*y) .* term;
    else
        G = 0;
    end
end

function G = Gprime_plus(y, z)
    % 计算 G'_+(y,z)，输入 y 和 z 必须为标量
    if y < z
        term = sqrt(z .* (1 - z) ./ (y .* (1 - y)));
        G = -1 + 2.*z + (1 - 2.*y) .* term;
    else
        G = 0;
    end
end

function G = G_minus(y, z)
    % 计算 G_-(y, z)，输入 y 和 z 必须为标量
    if y > 1 - z
        G = y + (1 - z).*(1 - 2.*y) - 2.*sqrt(z.*(1 - z).*y.*(1 - y));
    else
        G = 0;
    end
end

function G = G_plus(y, z)
    % 计算 G_+(y, z)，输入 y 和 z 必须为标量
    if y < z
        G = y + (1 - z).*(1 - 2.*y) + 2.*sqrt(z.*(1 - z).*y.*(1 - y));
    else
        G = 1;
    end
end

%%
function entropy = binary_entropy(p)
    if p < 0 || p > 1
        entropy = 0;
    end
    entropy = -p .* log2(p) - (1-p) .* log2(1-p);
end
