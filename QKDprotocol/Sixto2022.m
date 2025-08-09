%%
ncut = 30;
t = 4;
k = [0.5 0.2 10.^-4];
pk = [0.98,0.01,0.01];

% gamma_tilde_S = [0.5 0.51 0.503];
% gamma_tilde_D = [0.21 0.172 0.165];
% gamma_tilde_V = [10.^-4 10.^-4 10.^-4];
% 
% sigma_tilde_S = [0.032.*0.5 0.032.*0.51 0.034.*0.503];
% sigma_tilde_D = [0.07.*0.21 0.09.*0.172 0.091.*0.165];
% sigma_tilde_V = 10.^(-5).*gamma_tilde_V;

gamma_tilde_S = [0.5 0.5 0.5];
gamma_tilde_D = [0.2 0.2 0.2];
gamma_tilde_V = [10.^-4 10.^-4 10.^-4];

sigma_tilde_S = [0.032.*0.5 0.032.*0.5 0.034.*0.5];
sigma_tilde_D = [0.07.*0.2 0.09.*0.2 0.091.*0.2];
sigma_tilde_V = 10.^(-5).*gamma_tilde_V;

sigma_hat_S = sigma_tilde_S./gamma_tilde_S;
sigma_hat_D = sigma_tilde_D./gamma_tilde_D;
sigma_hat_V = sigma_tilde_V./gamma_tilde_V;

gamma_S = gamma_tilde_S;
gamma_D = gamma_tilde_D;
gamma_V = gamma_tilde_V;

sigma_S = invert_sigma(gamma_tilde_S, sigma_tilde_S, gamma_tilde_S - t.*sigma_tilde_S, gamma_tilde_S + t.*sigma_tilde_S);
sigma_D = invert_sigma(gamma_tilde_D, sigma_tilde_D, gamma_tilde_D - t.*sigma_tilde_D, gamma_tilde_D + t.*sigma_tilde_D);
sigma_V = invert_sigma(gamma_tilde_V, sigma_tilde_V, gamma_tilde_V - t.*sigma_tilde_V, gamma_tilde_V + t.*sigma_tilde_V);

gamma = [gamma_S; gamma_D; gamma_V];
sigma = [sigma_S'; sigma_D'; sigma_V'];
sigma_tilde = [sigma_tilde_S; sigma_tilde_D; sigma_tilde_V];

pd = 7.2.*10.^-8;
qZ = 0.99;
qX = 0.01;
L = 1:2:300;
alpha = 0.2;
eta_det = 0.65;

delta_A = 0.08;
f_EC = 1.16;
n = 0:ncut;
eta_ch = 10.^(-alpha.*L./10);
eta = eta_ch.*eta_det;

% 预先分配存储空间（如果需要在循环后使用结果）
infinite_key = zeros(length(L),1);

for i = 1:length(L)   
    % 计算部分
    [Dk_B, ek_B] = calculate_De_B(eta(i), pd, k, delta_A);
    y_tilde = zeros(1,length(n));
    h_tilde = zeros(1,length(n));
    for j = 1:length(n)
        [y_tilde(j), h_tilde(j)] = calculate_hy_tilde(delta_A, eta(i), pd, n(j));
    end
    pn = calculate_pn(gamma, sigma, sigma_tilde, t, n);
    CStao = caculate_CStao(pk, pn, n);
    [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
            = caculate_CScmts(CStao, y_tilde, h_tilde, n);
    
    % 线性规划部分
    [~, ZS1_bar_low] = linprog_sigma_y1h1(qZ, pk, pn, Dk_B, n, c_plus, c_minus, m_plus, m_minus);
    [~, XS1_bar_low] = linprog_sigma_y1h1(qX, pk, pn, Dk_B, n, c_plus, c_minus, m_plus, m_minus);
    [~, EXS1_bar_up] = linprog_sigma_h1h1(qX, pk, pn, ek_B, n, t_plus, t_minus, s_plus, s_minus);

    % 计算汇总值
    [Zkk_bar, ~] = caculate_ZXkk_bar(qZ, qX, pk, Dk_B);
    ZS_bar = sum(Zkk_bar(1,:));
    [EZkk_bar, EXkk_bar] = caculate_Ekk_bar(qZ, qX, pk, ek_B);
    EZS_bar = sum(EZkk_bar(1,:));

    % 存储结果
    infinite_key_temp = caculate_SKR(ZS1_bar_low, XS1_bar_low, EXS1_bar_up, ZS_bar, EZS_bar./ZS_bar, f_EC);
    if isempty(infinite_key_temp)
        infinite_key_temp = 0;
    end
    infinite_key(i) = infinite_key_temp;
end

%% 绘图
% Create figure with improved styling
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% Plot with logarithmic y-axis
semilogy(L, infinite_key, 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]);
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
function [nu_opt, fval] = linprog_sigma_y1h1(qZ, pk, pn, Dk_B, n, c_plus, c_minus, m_plus, m_minus)   
    % 高精度QKD线性规划问题求解
    
    % 目标函数构建 
    f = zeros(3*3*length(n), 1);
    for i = 1:3
        pos = (i-1)*length(n) + 2;  % n=1的位置
        f(pos) = qZ^2 * pk(1) * pk(i) * pn(1,i,2); 
    end
    
    % 使用更高精度的约束构建方式
    total_vars = 3*3*length(n);
    
    % 计算各约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    
    total_constr = num_constr1 + num_constr2 + num_constr3 + num_constr4;
    A = zeros(total_constr, total_vars);
    b = zeros(total_constr, 1);
    constr_idx = 1;
    
    % 约束1: Dk_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            base_pos = (a-1)*3*length(n) + (c-1)*length(n);
            for i = 1:length(n)
                row(base_pos + i) = pn(a,c,i);
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
            base_pos = (a-1)*3*length(n) + (c-1)*length(n);
            pn_sum = sum(pn(a,c,:));
            for i = 1:length(n)
                row(base_pos + i) = -pn(a,c,i);
            end
            A(constr_idx,:) = row;
            b(constr_idx) = 1 - pn_sum - Dk_B(a);
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
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n) + i;
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n) + i;
                    row(pos_a) = -m_plus(a,bb,c,i);
                    row(pos_bb) = 1;
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
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n) + i;
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n) + i;
                    row(pos_a) = m_minus(a,bb,c,i);
                    row(pos_bb) = -1;
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
function [nu_opt, fval] = linprog_sigma_h1h1(qX, pk, pn, ek_B, n, t_plus, t_minus, s_plus, s_minus)
    % 高精度QKD线性规划问题求解    
    
    % 目标函数构建 
    f = zeros(3*3*length(n), 1);
    for i = 1:3
        pos = (i-1)*length(n) + 2;  % n=1的位置
        f(pos) = qX^2 * pk(1) * pk(i) * pn(1,i,2); 
    end
    
    % 使用更高精度的约束构建方式
    total_vars = 3*3*length(n);    
    
    % 计算各约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    
    total_constr = num_constr1 + num_constr2 + num_constr3 + num_constr4;
    A = zeros(total_constr, total_vars);
    b = zeros(total_constr, 1);
    constr_idx = 1;
    
    % 约束1: ek_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            base_pos = (a-1)*3*length(n) + (c-1)*length(n);
            for i = 1:length(n)
                row(base_pos + i) = pn(a,c,i);
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
            base_pos = (a-1)*3*length(n) + (c-1)*length(n);
            pn_sum = sum(pn(a,c,:)); 
            for i = 1:length(n)
                row(base_pos + i) = -pn(a,c,i);
            end 
            A(constr_idx,:) = row;
            b(constr_idx) = 1 - pn_sum - ek_B(a);
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
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n) + i;
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n) + i;                    
                    row(pos_a) = -s_plus(a,bb,c,i);
                    row(pos_bb) = 1;                   
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
                    pos_a = (a-1)*3*length(n) + (c-1)*length(n) + i;
                    pos_bb = (bb-1)*3*length(n) + (c-1)*length(n) + i;                   
                    row(pos_a) = s_minus(a,bb,c,i);
                    row(pos_bb) = -1;                   
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
    %   n       : 光子数, 标量
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

function CStao = caculate_CStao(pk, pn, n) 
    % CStao 3x3x3
    CStao = zeros(3,3,3);
        for i = 1:3 % tao_a
            for j = 1:3 % tao_b
                for k =1:3  % tao_c
                    second_sigma = 0;
                        for r = 1:3
                            first_sigma = 0;
                                for t = 1:length(n)
                                    first_sigma = first_sigma + sqrt(pn(r,i,t).*pn(r,j,t));
                                end
                            second_sigma = second_sigma + pk(r).*first_sigma;
                        end
                    CStao(i,j,k) = second_sigma.^2;
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
function pn = calculate_pn(gamma, sigma, sigma_tilde, t, n)
    pn = zeros(3,3,length(n)); % pn是个3 x 3 x ncut+1的张量
    for i = 1:3
        for j = 1:3
            for k =1:length(n)
                pn(i,j,k) = local_compute(gamma(i,j), sigma(i,j), n(k), gamma(i,j) - t.*sigma_tilde(i,j), gamma(i,j) + t.*sigma_tilde(i,j));
            end
        end
    end
end

function p = local_compute(gamma, sigma, nk, lambda, Lambda)
    % 本地计算函数（处理标量参数）
    denominator = normal_cdf(gamma, sigma.^2, Lambda) - normal_cdf(gamma, sigma.^2, lambda);
    fact_ni = factorial(nk);
    integrand = @(a) normal_pdf(gamma, sigma.^2, a) .* exp(-a) .* a.^nk ./ (denominator.*fact_ni);
    p = integral(integrand, lambda, Lambda, 'RelTol', 1e-12);
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

function sigma = invert_sigma(gamma_tilde, sigma_tilde, lambda, Lambda)
    % 确保输入为列向量
    gamma_tilde = gamma_tilde(:);
    sigma_tilde = sigma_tilde(:);
    lambda = lambda(:);
    Lambda = Lambda(:);
    initial_guess = sigma_tilde .* 0.5;
    
    % 定义向量化方程（逐元素计算）
    fun = @(s) arrayfun(@(gt, st, l, L, s_val) ...
                        compute_sigma_tilde(gt, s_val, l, L) - st, ...
                        gamma_tilde, sigma_tilde, lambda, Lambda, s);
    
    % 使用 Levenberg-Marquardt 算法
    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
    sigma = abs(fsolve(fun, initial_guess, options));
end

function sigma_tilde = compute_sigma_tilde(gamma, sigma, lambda, Lambda)
    % 纯向量化计算
    alpha = (lambda - gamma) ./ sigma;
    beta = (Lambda - gamma) ./ sigma;
    
    phi_alpha = normal_pdf(0,1,alpha);
    phi_beta = normal_pdf(0,1,beta);
    Phi_alpha = normal_cdf(0,1,alpha);
    Phi_beta = normal_cdf(0,1,beta);
    
    term1 = (beta .* phi_beta - alpha .* phi_alpha) ./ (Phi_beta - Phi_alpha);
    term2 = ((phi_beta - phi_alpha) ./ (Phi_beta - Phi_alpha)).^2;
    
    sigma_tilde = abs(sqrt(sigma.^2 .* (1 - term1 - term2)));
end

%%
function pdf = normal_pdf(gamma, sigma2, x)
    % 计算正态分布的概率密度函数
    % gamma: 均值（原公式中的 y）
    % sigma2: 方差
    % x: 输入值
    sigma = sqrt(sigma2);
    pdf = (1 ./ (sigma .* sqrt(2 .* pi))) .* exp(-(x - gamma).^2 ./ (2 .* sigma2));
end

function cdf = normal_cdf(gamma, sigma2, x)
    % 计算正态分布的累积分布函数
    % gamma: 均值（原公式中的 y）
    % sigma2: 方差
    % x: 输入值
    sigma = sqrt(sigma2);
    z = (x - gamma) ./ sigma;  % 使用 ./ 按元素除法
    cdf = 0.5 .* (1 + erf(z ./ sqrt(2)));
end

%%
function entropy = binary_entropy(p)
%     if p < 0 || p > -1
%         entropy = 0;
%     end
    entropy = -p .* log2(p) - (1-p) .* log2(1-p);
end
