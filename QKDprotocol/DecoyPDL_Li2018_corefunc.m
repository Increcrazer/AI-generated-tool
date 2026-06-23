% 参考：
%   C. Li et al., Phys. Rev. A 98, 042324 (2018), Fig. 4(a).

clear; clc;

%% Li2018 Table I 参数
ed = 0.015;
Y0 = 2e-5;
eta_Bob = 0.2;
alpha_dB_km = 0.2;
f_EC = 1.16;
e0 = 0.5;
q = 0.25;                   % 按总发射脉冲归一化：四态等概率 + 标准 BB84 筛基

PDL_list_dB = [0, 1.6, 3, 5, 10];
distance_km = 0:1:150;

% 每个距离、每个 PDL 单独优化 mu。Li2018 中大 PDL 时 mu_opt 可超过 1。
mu_bounds = [1e-4, 5];

R_post = zeros(numel(PDL_list_dB),numel(distance_km));
R_no_post = zeros(numel(PDL_list_dB),numel(distance_km));
mu_opt_post = zeros(size(R_post));
P_opt_post = zeros(size(R_post));
mu_opt_no_post = zeros(size(R_post));

%% 逐 PDL 和距离计算 Fig. 4(a)
for ip = 1:numel(PDL_list_dB)
    PDL_dB = PDL_list_dB(ip);
    L_pdl = 10^(-PDL_dB/10);

    for il = 1:numel(distance_km)
        eta_sys = eta_Bob * 10^(-alpha_dB_km*distance_km(il)/10);

        % Fig. 4(a): P = P_opt(mu)，并对 mu 优化。
        obj_post = @(mu) -li2018_rate(mu,L_pdl,eta_sys,Y0,ed,e0,f_EC,q,true);
        [mu_star_post,neg_rate_post] = fminbnd(obj_post,mu_bounds(1),mu_bounds(2), ...
            optimset('TolX',1e-8,'Display','off'));
        R_post(ip,il) = max(0,-neg_rate_post);
        mu_opt_post(ip,il) = mu_star_post;
        P_opt_post(ip,il) = li2018_p_opt(mu_star_post,L_pdl);

        % Li2018 Fig. 4(b) / Table II 的对照：保持同一个 mu_s,opt，只令 P=1。
        R_no_post(ip,il) = max(0,li2018_rate(mu_star_post,L_pdl,eta_sys,Y0,ed,e0,f_EC,q,false));
        mu_opt_no_post(ip,il) = mu_star_post;
    end
end

%% 画 Li2018 Fig. 4(a) 对应的后选择曲线
figure('Color','white','Position',[100 100 850 620]);
hold on;
colors = lines(numel(PDL_list_dB));
for ip = 1:numel(PDL_list_dB)
    semilogy(distance_km,max(R_post(ip,:),1e-15),'LineWidth',2.0, ...
        'Color',colors(ip,:),'DisplayName',sprintf('PDL = %.1f dB',PDL_list_dB(ip)));
end
grid on;
set(gca,'YScale','log');
xlabel('Distance (km)');
ylabel('Secret key rate (bit/pulse)');
title('Li2018 Fig. 4(a) check: four-state BB84 with postselection');
ylim([1e-8 1]);
xlim([0 max(distance_km)]);
legend('Location','southwest');

%% 额外诊断图：后选择 vs 无后选择
figure('Color','white','Position',[160 120 850 620]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
for ip = 1:numel(PDL_list_dB)
    semilogy(distance_km,max(R_no_post(ip,:),1e-15),'LineWidth',2.0, ...
        'Color',colors(ip,:),'DisplayName',sprintf('PDL = %.1f dB',PDL_list_dB(ip)));
end
grid on;
set(gca,'YScale','log');
xlabel('Distance (km)');
ylabel('Secret key rate (bit/pulse)');
title('No postselection, P = 1');
ylim([1e-8 1]);
xlim([0 max(distance_km)]);

nexttile;
hold on;
for ip = 1:numel(PDL_list_dB)
    plot(distance_km,P_opt_post(ip,:),'LineWidth',2.0, ...
        'Color',colors(ip,:),'DisplayName',sprintf('PDL = %.1f dB',PDL_list_dB(ip)));
end
grid on;
xlabel('Distance (km)');
ylabel('P_{opt}');
title('Postselection probability');
ylim([0 1.05]);
legend('Location','southwest');

%% 打印 80 km 处的数值，和 Li2018 Table II 做量级/趋势对照
target_distance = 80;
[~,idx80] = min(abs(distance_km-target_distance));
fprintf('Li2018 asymptotic four-state PDL check at %.0f km\n',distance_km(idx80));
fprintf('PDL(dB)    R_post        R_no_post     increase       mu_post    P_opt\n');
for ip = 1:numel(PDL_list_dB)
    inc = R_post(ip,idx80)./max(R_no_post(ip,idx80),realmin)-1;
    fprintf('%6.1f   %.4e   %.4e   %10.3g   %.4f   %.4f\n', ...
        PDL_list_dB(ip),R_post(ip,idx80),R_no_post(ip,idx80), ...
        inc,mu_opt_post(ip,idx80),P_opt_post(ip,idx80));
end

%% 局部函数
function R = li2018_rate(mu,L_pdl,eta_sys,Y0,ed,e0,f_EC,q,use_postselection)
%LI2018_RATE Li2018 无限诱骗渐近密钥率。
% mu 是未受 PDL 损耗的 V/A 信号强度；H/D 信号强度为 L_pdl*mu。
% use_postselection=true 时使用 P=P_opt；false 时 P=1。
    mu_H = L_pdl*mu;
    mu_V = mu;

    if use_postselection
        P = li2018_p_opt(mu,L_pdl);
    else
        P = 1;
    end

    Y1 = 1 - (1-Y0)*(1-eta_sys);
    e1_phase = (e0*Y0 + ed*(Y1-Y0))/Y1;  % 四态 BB84: 用 X 基单光子误码估计相位误码

    Q_H = gain_coherent(mu_H,eta_sys,Y0);
    Q_V = gain_coherent(mu_V,eta_sys,Y0);
    EQ_H = error_gain_coherent(mu_H,eta_sys,Y0,ed,e0);
    EQ_V = error_gain_coherent(mu_V,eta_sys,Y0,ed,e0);

    Qs_tilde = 0.5*Q_H + 0.5*P*Q_V;
    Es_tilde = (0.5*EQ_H + 0.5*P*EQ_V)/max(Qs_tilde,realmin);

    Q1_tilde = min(mu_H*exp(-mu_H), P*mu_V*exp(-mu_V))*Y1;

    R = q*(Q1_tilde*(1-binary_entropy(e1_phase)) ...
        - Qs_tilde*f_EC*binary_entropy(Es_tilde));
end

function P = li2018_p_opt(mu,L_pdl)
%LI2018_P_OPT Li2018 Eq. (18): L*mu*exp(-L*mu) = P*mu*exp(-mu)。
% 因此 P = L*exp((1-L)*mu)，并限制在 [0,1] 内。
    P = L_pdl*exp((1-L_pdl)*mu);
    P = min(1,max(0,P));
end

function Q = gain_coherent(mu,eta_sys,Y0)
%GAIN_COHERENT 相干态总点击增益 Q_mu。
    Q = 1 - (1-Y0)*exp(-eta_sys*mu);
end

function EQ = error_gain_coherent(mu,eta_sys,Y0,ed,e0)
%ERROR_GAIN_COHERENT 相干态错误增益 E_mu Q_mu。
    EQ = e0*Y0 + ed*(1-exp(-eta_sys*mu))*(1-Y0);
end

function h = binary_entropy(x)
%BINARY_ENTROPY 二元 Shannon 熵。
    x = min(max(x,realmin),1-realmin);
    h = -x*log2(x) - (1-x)*log2(1-x);
end
