Pout_dBm = -26; % 输入功率 (dBm)
Pout = 10^(Pout_dBm/10);
lambda = 1550 * 10^(-9); % 波长 (m)
sp_energy = 6.626e-34 * 3e8 / lambda; % 单光子能量 (J)
rou_1530 = 3.2e-9; % 拉曼散射系数 (W⁻¹·km⁻¹)
eta = 0.7; % 探测效率
alpha_lg = 0.2; % 光纤损耗 (dB/km)
alpha_ln = log(10)/10 * alpha_lg; % 转换为线性单位 (Np/km)
L = 0:100; % 光纤长度（km）

% 计算前向和后向拉曼散射功率 (mW)
Pram_f_power = Pout .* L * rou_1530 * eta; 
Pram_b_power = Pout .* sinh(alpha_ln * L) / alpha_ln * rou_1530 * eta;

% 计算光子数率 (光子数/秒)
Pram_f_photon = Pram_f_power * 1e-3 / sp_energy; % 功率 (mW → W) → 光子数/秒
Pram_b_photon = Pram_b_power * 1e-3 / sp_energy;

% 创建双纵坐标轴图形
figure;
yyaxis left; % 左侧轴：功率 (mW)
semilogy(L, Pram_f_power, 'LineWidth', 2, 'Color', [0, 0.5, 0.8], 'DisplayName', 'Forward Raman (Power)');
hold on;
semilogy(L, Pram_b_power, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'Backward Raman (Power)');
ylabel('Raman Power (mW)', 'FontSize', 12);

yyaxis right; % 右侧轴：光子数率 (光子数/秒)
plot(L, Pram_f_photon, '--', 'LineWidth', 2, 'Color', [0, 0.7, 0.3], 'DisplayName', 'Forward Raman (Photon Rate)');
plot(L, Pram_b_photon, '--', 'LineWidth', 2, 'Color', [1, 0.5, 0], 'DisplayName', 'Backward Raman (Photon Rate)');
ylabel('Photon Rate (photons/s)', 'FontSize', 12);
hold off;

% 美化图像
title('Raman Scattering', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Fiber Length (km)', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
box on;

% 调整坐标轴范围
xlim([1, 100]);
yyaxis left;
ylim([0, max([Pram_f_power, Pram_b_power]) * 1.1]); % 左侧Y轴范围
yyaxis right;
ylim([0, max([Pram_f_photon, Pram_b_photon]) * 1.1]); % 右侧Y轴范围

% 设置坐标轴字体大小
set(gca, 'FontSize', 10);
