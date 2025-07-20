  %%
    % 输入：
    % epsilon_cor, epsilon_sec: security criteria
    % qX: 选X基概率(efficient BB84), N: A发送的块长, f[Hz]: 频率, gate_width[s]: 门宽,
    % 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[m]: A到B距离, alpha[dB/m]: 光路损耗
    % ita_Bob_detect: Bob解码端总效率（假设通过PC使得探测器效率几乎一致）,包括滤波解码端超导等，不包括门宽
    % e_mis_X(Z): A端态制备的X(Z)基错误率，用X(Z)基对比度反算得到, f_EC: 纠错效率
    % p_ap: 后脉冲概率, dc_count[/s]: 超导单通道暗计数, 是个1x4向量（四个超导）,deadtime[s]: 死时间
    % 输出：
    % R_bitperpulse, R_bitpersecond:密钥率
    % e_obs: X基比特误码率（用于成码的是X基）
    % phi_X: 估计得到的用于成码的X基的相位误码率
    % nX: B接收到的X基块长
    
%% 输入参数
epsilon_cor = 10^-15;
epsilon_sec = 10^-15;
qX = 0.5;
N = 10^20;
f = 1.25*10^5; 
gate_width = 1/f; 
width_3dB = 70*10^(-12);
k = [0.6 0.1 0.0002]; 
pk = [0.5 0.25 0.25];
L = 0:2:200; 
alpha = 0.2; 
ita_Bob_detect = 0.1;
e_mis_X = 0.005;
e_mis_Z = 0.005;
f_EC = 1.16; 
p_ap = 0.04;
dc_count = [f*6*10^(-7) f*6*10^(-7) f*6*10^(-7) f*6*10^(-7)];
deadtime = 0;

%%
R_bitperpulse = zeros(1,length(L));
e_obs = zeros(1,length(L));
phi_X = zeros(1,length(L));
for i = 1:length(L)
    [R_bitperpulse(i), R_bitpersecond, e_obs(i), phi_X(i), nX, nuZ_1] = Decoy_Lim2014_corefunc ... 
                (epsilon_cor, epsilon_sec, ...
                 qX, N, f, gate_width, width_3dB, ...
                 k, pk, ...
                 L(i), alpha, ita_Bob_detect, ...
                 e_mis_X, e_mis_Z, f_EC, ...
                 p_ap, dc_count, deadtime);
    if R_bitperpulse(i)<0
        R_bitperpulse(i)=0;
    end
end

%% 绘图
% Create figure with improved styling
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% Plot with logarithmic y-axis
semilogy(L, R_bitperpulse, 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]);
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
