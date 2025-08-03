  %%
    % 输入：
    % qZ: 选Z基概率(efficient BB84), ncut: 诱骗态截止光子数, gate_width[s]: 门宽, 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[m]: A到B距离, alpha[dB/m]: 光路损耗
    % eta_Bob_detect: Bob解码端总效率（假设通过PC使得探测器效率几乎一致）,包括滤波解码端超导等，不包括门宽
    % delta_A: 偏振基偏差角度, f_EC: 纠错效率
    % pd: 暗计数概率, delta_max: 最大偏差
    % 输出：
    % R_bitperpulse:密钥率
    
%% 输入参数
qZ = 0.95;
ncut = 20;
f = 1.25*10^5; 
gate_width = 1/f; 
width_3dB = 70*10^(-12);
k = [0.5 0.2 10.^-4];
pk = [0.95 0.03 0.02];
L = 0:1:150; 
alpha = 0.2; 
eta_Bob_detect = 0.65;
delta_A = 0.08;
f_EC = 1.16; 
pd = 7.2.*10.^-8;
delta_max = [10^-2 10^-3 10^-4];

%%
R_bitperpulse = zeros(length(delta_max),length(L));
for i = 1:length(delta_max)
    for j = 1:length(L)
        R_bitperpulse(i,j) = Decoy_Sixto2022MI_corefunc ... 
                (qZ, ncut, gate_width, width_3dB, ...
                 k, pk, ...
                 L(j), alpha, eta_Bob_detect, ...
                 delta_A, f_EC, ...
                 pd, delta_max(i));
        if R_bitperpulse(j)<0
            R_bitperpulse(j)=0;
        end
    end
end

%% 绘图
% Create figure with improved styling
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% Plot with logarithmic y-axis
semilogy(L, R_bitperpulse(1,:), 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'delta max = 1E-2'); hold on;
semilogy(L, R_bitperpulse(2,:), 'LineWidth', 2.5, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'delta max = 1E-3'); 
semilogy(L, R_bitperpulse(3,:), 'LineWidth', 2.5, 'Color', [0.9290, 0.6940, 0.1250], 'DisplayName', 'delta max = 1E-4'); 
grid on;

% Add labels and title with larger fonts
xlabel('Distance (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Secret Key Rate (bits/pulse)', 'FontSize', 14, 'FontWeight', 'bold');
title('QKD Key Rate vs Distance', 'FontSize', 16, 'FontWeight', 'bold');

% Add legend
legend('show', 'FontSize', 12, 'Location', 'best');

% Adjust axes properties
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.YScale = 'log';
ax.YMinorTick = 'on';
ax.XMinorTick = 'on';
ax.GridAlpha = 0.3;
