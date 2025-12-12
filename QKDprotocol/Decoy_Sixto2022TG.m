%% 并非原始Sixto2022 TG，而是使用Trefilov2025的L+1距离使用L距离的输出yn进行迭代的思想，由于距离之间存在参数传递，所以无法写成corefunc形式
%% Sixto2022 TG的FIG.3由于未知原因复现不出，大概率是作者代码写错了
%% 该版本经师弟yjk修改跑通后进行修订整理

%% user parameter set
% gamma_tilde_S = [0.5 0.51 0.503];
% gamma_tilde_D = [0.21 0.172 0.165];
% gamma_tilde_V = [10.^-4 10.^-4 10.^-4];
% 
% sigma_tilde_S = [0.032.*0.5 0.032.*0.51 0.034.*0.503];
% sigma_tilde_D = [0.07.*0.21 0.09.*0.172 0.091.*0.165];
% sigma_tilde_V = 10.^(-5).*gamma_tilde_V;

%% literature parameter set, to test the validity of the code
%% Sixto2022 
%% 有相关性
% gamma_tilde_S = [0.5 0.51 0.503];
% gamma_tilde_D = [0.21 0.172 0.165];
% gamma_tilde_V = [10.^-4 10.^-4 10.^-4];
% 
% sigma_tilde_S = [0.032.*0.5 0.032.*0.51 0.034.*0.503];
% sigma_tilde_D = [0.07.*0.21 0.09.*0.172 0.091.*0.165];
% sigma_tilde_V = 10.^(-5).*gamma_tilde_V;
%% 无相关性
% gamma_tilde_S = [0.5 0.5 0.5];
% gamma_tilde_D = [0.2 0.2 0.2];
% gamma_tilde_V = [10.^-4 10.^-4 10.^-4];
% 
% sigma_tilde_S = [0.032.*0.5 0.032.*0.5 0.034.*0.5];
% sigma_tilde_D = [0.07.*0.2 0.09.*0.2 0.091.*0.2];
% sigma_tilde_V = 1e-5.*gamma_tilde_V;

%% Trefilov2025
%% System A（无相关性）
% gamma_tilde_S = [0.638635 0.638635 0.638635];
% gamma_tilde_D = [0.276372 0.276372 0.276372];
% gamma_tilde_V = [0.042044 0.042044 0.042044];
% 
% sigma_tilde_S = [0.025245 0.025245 0.025245];
% sigma_tilde_D = [0.011834 0.011834 0.011834];
% sigma_tilde_V = [0.011438 0.011438 0.011438];

%% System A（有相关性）
gamma_tilde_S = [0.639209 0.639251 0.635059];
gamma_tilde_D = [0.278423 0.267539 0.273149];
gamma_tilde_V = [0.042104 0.042053 0.041681];

sigma_tilde_S = [0.025230 0.025164 0.024665];
sigma_tilde_D = [0.011295 0.011064 0.010465];
sigma_tilde_V = [0.011465 0.011400 0.011308];

%% 2025 System B
%% System B 归一化（无相关）
% gamma_tilde_S0 = [1 1 1];
% gamma_tilde_D0 = [0.331205 0.331205 0.331205];
% gamma_tilde_V0 = [0.080408 0.080408 0.080408];
% 
% sigma_tilde_S0 = [0.03082 0.03082 0.03082];
% sigma_tilde_D0 = [0.031052 0.031052 0.031052];
% sigma_tilde_V0 = [0.019021 0.019021 0.019021];
%% System B 归一化（有相关）
% gamma_tilde_S0 = [0.996082 1.004939 1.002882];
% gamma_tilde_D0 = [0.32323 0.339161 0.339283];
% gamma_tilde_V0 = [0.081123 0.080498 0.078888];
% 
% sigma_tilde_S0 = [0.030962 0.030197 0.030115];
% sigma_tilde_D0 = [0.030018 0.030004 0.029979];
% sigma_tilde_V0 = [0.01879 0.019254 0.019157];
% %% 真实值
% % mu = 0.43/ 0.22/ 0.11
% mu = 0.11;
% gamma_tilde_S = mu.*gamma_tilde_S0;
% gamma_tilde_D = mu.*gamma_tilde_D0;
% gamma_tilde_V = mu.*gamma_tilde_V0;
% sigma_tilde_S = mu.*sigma_tilde_S0;
% sigma_tilde_D = mu.*sigma_tilde_D0;
% sigma_tilde_V = mu.*sigma_tilde_V0;

%% 生成距离（线性/对数）
nL = 100;
l = linspace(0,1,nL);
Lmax = 180;%信道长度
% L = linspace(0,Lmax,nL);
dens_var = 50;%数据点稀疏程度控制
L = Lmax * log10(1 + dens_var*l)/log10(1+dens_var);

%% 参数初始化
ncut = 20;
t = 4;
pk = [0.5,0.3,0.2];

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
alpha = 0.2;    
eta_det = 0.65; 

delta_A = 0.08; 
f_EC = 1.16;    
n = 0:ncut;
eta_ch = 10.^(-alpha.*L./10);  
eta = eta_ch.*eta_det;  

% 预先分配存储空
infinite_key = zeros(length(L),1);

%% 计算部分
for i = 1:length(L)   
    [Dk_B, ek_B] = calculate_De_B(eta(i), pd, gamma, sigma, sigma_tilde, t, delta_A);  

    [y_tilde, h_tilde] = calculate_hy_tilde(delta_A, eta(i), pd, n);%计算参考产出率和误码率

    pn = calculate_pn(gamma, sigma, sigma_tilde, t, n);
    CStao = caculate_CStao(pk, pn, n);
    
    if i==1
    [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
            = caculate_CScmts_ini(CStao, y_tilde, h_tilde, n);
    else
    [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
            = caculate_CScmts_iterate(CStao, yn_ZS, hn_EXS, n);
    end
    % 线性规划部分
    [yn_ZS, ZS1_bar_low, ~] = linprog_sigma_y1h1(qZ, pk, pn, Dk_B, n, c_plus, c_minus, m_plus, m_minus);
    XS1_bar_low = ZS1_bar_low*qX^2/qZ^2;
    [hn_EXS, EXS1_bar_up, ~] = linprog_sigma_h1h1(qX, pk, pn, ek_B, n, t_plus, t_minus, s_plus, s_minus);%求解线性规划得到单光子错误率的上界
    
    % 对y和h进行修正，使之不为0和1
    yn_ZS(yn_ZS>=1) = 1-1e-16;
    yn_ZS(yn_ZS<=0) = 1e-16;
    hn_EXS(hn_EXS>=1) = 1-1e-16;
    hn_EXS(hn_EXS<=0) = 1e-16;

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
%Create figure with improved styling
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% Plot with logarithmic y-axis
semilogy(L, infinite_key, 'LineWidth', 2.5);
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

%% infinite_key
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
function [nu_opt, fval, exitflag] = linprog_sigma_y1h1(qZ, pk, pn, Dk_B, n, c_plus, c_minus, m_plus, m_minus)   
    % 总变量数量
    total_vars = 3*3*length(n);    

    % 计算各约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    num_constr5 = 3*3*(length(n)-1);%单调性约束
    total_constr = num_constr1 + num_constr2 + num_constr3 + num_constr4 + num_constr5;

    % 目标函数构建 
    f = zeros(total_vars, 1);
    for c = 1:3
        pos = pos_convert(1, c, 1+1, length(n));  % n=1的位置
        f(pos) = qZ^2 * pk(1) * pk(c) * pn(1,c,2); 
    end
    
    A = zeros(total_constr, total_vars);
    b = zeros(total_constr, 1);
    constr_idx = 1;

    % 约束1: Dk_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            for i = 1:length(n)
                pos = pos_convert(a,c,i, length(n));
                row(pos) = pn(a,c,i);
            end           
            A(constr_idx,:) = row;
            b(constr_idx) = Dk_B(a,c);
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束2: 1-Dk_B约束 
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pn_sum = sum(pn(a,c,:));
            for i = 1:length(n)
                pos = pos_convert(a,c,i, length(n));
                row(pos) = -pn(a,c,i);
            end
            A(constr_idx,:) = row;
            b(constr_idx) = 1 - pn_sum - Dk_B(a,c);
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
                    pos_a = pos_convert(a,c,i, length(n));
                    pos_bb = pos_convert(bb,c,i, length(n));
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
                    pos_a = pos_convert(a,c,i, length(n));
                    pos_bb = pos_convert(bb,c,i, length(n));
                    row(pos_a) = m_minus(a,bb,c,i);
                    row(pos_bb) = -1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = -c_minus(a,bb,c,i);
                    constr_idx = constr_idx + 1;
                end
            end
        end
    end
    
    % 约束5: 单调性约束 y(a,c,k+1)>=y(a,c,k)
    for a = 1:3
        for c = 1:3
            for i = 1:length(n)-1
                    row = zeros(1, total_vars);
                    pos_k = pos_convert(a,c,i, length(n));
                    pos_kplus = pos_convert(a,c,i+1, length(n));
                    row(pos_k) = 1;
                    row(pos_kplus) = -1;
                    A(constr_idx,:) = row;
                    b(constr_idx) = 0;
                    constr_idx = constr_idx + 1;
            end
        end
    end
    
    % 高精度求解选项
    options = optimoptions('linprog', ...
        'Algorithm', 'dual-simplex', ...  % 保持对偶单纯形法
        'ConstraintTolerance', 1e-9, ...  % 约束容差
        'OptimalityTolerance', 1e-10, ...  % 最优性容差
        'MaxIterations', 100000, ...      % 最大迭代次数提高到10万次
        'Display', 'final');               % 显示迭代过程便于调试
    
    % 求解线性规划
    [nu_opt, fval, exitflag] = linprog(f, A, b, [], [], zeros(total_vars,1), ones(total_vars,1), options);
end

function [nu_opt, fval, exitflag] = linprog_sigma_h1h1(qX, pk, pn, ek_B, n, t_plus, t_minus, s_plus, s_minus)
    % 总变量数量
    total_vars = 3*3*length(n);    
    
    % 计算各约束数量
    num_constr1 = 3*3;  % 约束1的数量
    num_constr2 = 3*3;  % 约束2的数量
    num_constr3 = 3*3*2*length(n); % 约束3的数量 (2因为bb≠a)
    num_constr4 = 3*3*2*length(n); % 约束4的数量
    total_constr = num_constr1 + num_constr2 + num_constr3 + num_constr4;

    % 目标函数构建 
    f = zeros(3*3*length(n), 1);
    for c = 1:3
        pos = pos_convert(1,c,1+1,length(n));  % n=1的位置
        f(pos) = qX^2 * pk(1) * pk(c) * pn(1,c,1+1); 
    end
    
    A = zeros(total_constr, total_vars);
    b = zeros(total_constr, 1);
    constr_idx = 1;
    
    % 约束1: ek_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            for i = 1:length(n)
                pos = pos_convert(a,c,i,length(n));
                row(pos) = pn(a,c,i);
            end
            A(constr_idx,:) = row;
            b(constr_idx) = ek_B(a,c);
            constr_idx = constr_idx + 1;
        end
    end
    
    % 约束2: 1-ek_B约束
    for a = 1:3
        for c = 1:3
            row = zeros(1, total_vars);
            pn_sum = sum(pn(a,c,:)); 
            for i = 1:length(n)
                pos = pos_convert(a,c,i,length(n));
                row(pos) = -pn(a,c,i);
            end 
            A(constr_idx,:) = row;
            b(constr_idx) = 1 - pn_sum - ek_B(a,c);
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
                    pos_a = pos_convert(a,c,i,length(n));
                    pos_bb = pos_convert(bb,c,i,length(n));             
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
                    pos_a = pos_convert(a,c,i,length(n));
                    pos_bb = pos_convert(bb,c,i,length(n));                     
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
        'ConstraintTolerance', 1e-9, ...  % 约束容差
        'OptimalityTolerance', 1e-10, ...  % 最优性容差
        'MaxIterations', 100000, ...      % 最大迭代次数
        'Display', 'off');              
    
    % 求解线性规划
    [nu_opt, fval, exitflag] = linprog(-f, A, b, [], [], zeros(total_vars,1), ones(total_vars,1),  options);
    fval = -fval;
end

%%
function pos = pos_convert(a, c, i, len_n)
    pos = (a-1).*3.*len_n + (c-1).*len_n + i;
end

%%
function [Dk_B, ek_B] = calculate_De_B(eta, pd, gamma, sigma, sigma_tilde, t, delta_A)
    % 计算探测率Dk_B和错误率ek_B（支持k为向量输入）
    % 输入参数：
    %   eta     : 系统总效率 (η = η_def .* η_ch)
    %   pd      : 探测器暗计数概率
    %   gamma   : 每种模式的平均强度向量 3x3
    %   sigma   : 每种模式的强度标准差向量 3x3
    %   delta_A : 偏振失配角 (rad)
    % 输出：
    %   Dk_B    : B端探测率向量 3x3
    %   ek_B    : B端错误率向量 3x3

    % 计算探测率 Dk_B (公式类似文中y_tilde)
     Dk_B = zeros(3,3);
     ek_B = zeros(3,3);
     for i=1:3
         for j=1:3
             Dk_B(i,j) = local_compute1(gamma(i,j), sigma(i,j), gamma(i,j) - t.*sigma_tilde(i,j), gamma(i,j) + t.*sigma_tilde(i,j), eta, pd);
             ek_B(i,j) = local_compute2(gamma(i,j), sigma(i,j), gamma(i,j) - t.*sigma_tilde(i,j), gamma(i,j) + t.*sigma_tilde(i,j), eta, pd, delta_A);
         end
     end
 end

function Dk_B = local_compute1(gamma, sigma, lambda, Lambda, eta, pd)
    % 本地计算函数（处理标量参数）
    denominator = normal_cdf(gamma, sigma.^2, Lambda) - normal_cdf(gamma, sigma.^2, lambda);
    g = @(a) normal_pdf(gamma, sigma.^2, a)  ./ (denominator);
    f = @(a) g(a).*(1 - (1 - pd).^2 .* exp(-eta .* a));
    Dk_B = integral(f, lambda, Lambda, 'RelTol', 1e-12,'AbsTol', 1e-15);
end

function ek_B = local_compute2(gamma, sigma, lambda, Lambda, eta, pd, delta_A)
    % 本地计算函数（处理标量参数）
    denominator = normal_cdf(gamma, sigma.^2, Lambda) - normal_cdf(gamma, sigma.^2, lambda);
    g = @(a) normal_pdf(gamma, sigma.^2, a)  ./ (denominator);
    h = @(a) (exp(-eta .* a .* cos(delta_A).^2) - exp(-eta .* a .* sin(delta_A).^2)) ./ 2;
    f = @(a) g(a).*(pd.^2 ./ 2 + pd .* (1 - pd) .* (1 + h(a)) + (1 - pd).^2 .* (0.5 + h(a) - 0.5 .* exp(-eta .* a)));
    ek_B = integral(f, lambda, Lambda, 'RelTol', 1e-12,'AbsTol', 1e-15);
end

%%
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

%%
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
                    % 相当重要的修正
                    if (second_sigma > 1)
                        CStao(i,j,k) = 1;
                    else
                        CStao(i,j,k) = second_sigma.^2;
                    end
                end
            end
        end
end

%% 
function [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
        = caculate_CScmts_ini(CStao, y_tilde, h_tilde, n)
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

function [c_plus, c_minus, m_plus, m_minus, t_plus, t_minus, s_plus, s_minus] ...
        = caculate_CScmts_iterate(CStao, y_tilde, h_tilde, n)
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
                    c_plus(i,j,k,l) = G_plus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)) - Gprime_plus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)).*y_tilde(pos_convert(i,k,l,length(n)));
                    c_minus(i,j,k,l) = G_minus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)) - Gprime_minus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)).*y_tilde(pos_convert(i,k,l,length(n)));
                    m_plus(i,j,k,l) = Gprime_plus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k));
                    m_minus(i,j,k,l) = Gprime_minus(y_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k));
                    t_plus(i,j,k,l) = G_plus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)) - Gprime_plus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)).*h_tilde(pos_convert(i,k,l,length(n)));
                    t_minus(i,j,k,l) = G_minus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)) - Gprime_minus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k)).*h_tilde(pos_convert(i,k,l,length(n)));
                    s_plus(i,j,k,l) = Gprime_plus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k));
                    s_minus(i,j,k,l) = Gprime_minus(h_tilde(pos_convert(i,k,l,length(n))),CStao(i,j,k));                    
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

%% 
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
