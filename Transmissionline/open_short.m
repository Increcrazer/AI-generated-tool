% 传输线电压和电流波动画 - 短路和开路情况
clear; close all; clc;

% 参数设置
Z0 = 50;            % 传输线特性阻抗 (Ω)
f = 1e9;            % 频率 1GHz
c = 3e8;            % 光速 (m/s)
lambda = c/f;       % 波长 (m)
beta = 2*pi/lambda; % 相位常数 (rad/m)
L = 2*lambda;       % 传输线长度 (m)
z = linspace(0, L, 1000); % 传输线位置坐标
t = linspace(0, 3/f, 100); % 时间点 (3个周期)

% 创建图形窗口
figure('Position', [100, 100, 1200, 800]);
sgtitle('传输线电压和电流波动画 (Z_0 = 50Ω)');

% 情况1: 短路负载 (ZL = 0)
ZL_short = 0;
[V_short, I_short] = transmissionLineWave(Z0, ZL_short, beta, z);

% 情况2: 开路负载 (ZL -> ∞)
ZL_open = 1e10;
[V_open, I_open] = transmissionLineWave(Z0, ZL_open, beta, z);

% 创建动画
for k = 1:length(t)
    % 计算瞬时值
    V_short_t = real(V_short * exp(1j*2*pi*f*t(k)));
    I_short_t = real(I_short * exp(1j*2*pi*f*t(k)));
    V_open_t = real(V_open * exp(1j*2*pi*f*t(k)));
    I_open_t = real(I_open * exp(1j*2*pi*f*t(k)));
    
    % 短路情况
    subplot(2,2,1);
    plot(z/lambda, V_short_t, 'b', 'LineWidth', 2);
    ylim([-2.2 2.2]);
    title('短路负载 - 电压波');
    xlabel('位置 z/\lambda'); ylabel('电压 (V)');
    grid on;
    
    subplot(2,2,2);
    plot(z/lambda, I_short_t, 'r', 'LineWidth', 2);
    ylim([-0.05 0.05]);
    title('短路负载 - 电流波');
    xlabel('位置 z/\lambda'); ylabel('电流 (A)');
    grid on;
    
    % 开路情况
    subplot(2,2,3);
    plot(z/lambda, V_open_t, 'b', 'LineWidth', 2);
    ylim([-2.2 2.2]);
    title('开路负载 - 电压波');
    xlabel('位置 z/\lambda'); ylabel('电压 (V)');
    grid on;
    
    subplot(2,2,4);
    plot(z/lambda, I_open_t, 'r', 'LineWidth', 2);
    ylim([-0.05 0.05]);
    title('开路负载 - 电流波');
    xlabel('位置 z/\lambda'); ylabel('电流 (A)');
    grid on;
    
    % 添加时间标记
    annotation('textbox', [0.15, 0.9, 0.1, 0.1], 'String',...
        sprintf('时间: %.2f ns', t(k)*1e9), 'EdgeColor', 'none');
    
    drawnow;
    pause(0.05);
    
    frame = getframe(gcf);
end

% 传输线波计算函数

function [V, I] = transmissionLineWave(Z0, ZL, beta, z)
    % 计算反射系数
    Gamma = (ZL - Z0)/(ZL + Z0);
    
    % 假设入射波幅度为1
    V_inc = 1;
    I_inc = V_inc / Z0;
    
    % 计算总电压和电流分布
    V = V_inc * (exp(-1j*beta*z) + Gamma * exp(1j*beta*z));
    I = I_inc * (exp(-1j*beta*z) - Gamma * exp(1j*beta*z));
end
