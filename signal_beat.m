%% 功能：两路光信号的拍频效应及理论相位分析（可选相位突变）
clear; clc;

%% 参数设置
fs = 1e6;           % 采样频率 (Hz)
T = 1e-3;           % 信号持续时间 (s)
t = 0:1/fs:T-1/fs;  % 时间向量

f1 = 100e3;         % 信号1频率 (Hz)
f2 = 105e3;         % 信号2频率 (Hz)
A1 = 1;             % 信号1幅度
A2 = 1;             % 信号2幅度
phi1 = 0;           % 信号1初始相位
phi2 = 0;           % 信号2初始相位（默认无突变）

%% 可选：在t=0.5ms时对signal2添加0.3π相位突变（取消注释启用）
% phase_step = 0.3 * pi;
% step_idx = find(t >= 0.5e-3, 1);
% phi2 = zeros(size(t));
% phi2(step_idx:end) = phase_step;

%% 生成信号
signal1 = A1 * sin(2*pi*f1*t + phi1);
signal2 = A2 * sin(2*pi*f2*t + phi2);
combined_signal = signal1 + signal2;

%% 理论拍频相位计算
% 公式：φ_beat = (φ1 - φ2)/2 + π(f1 - f2)t
beat_phase = (phi1 - phi2)/2 + pi*(f1 - f2)*t;

%% 绘图
figure;

% 子图1：原始信号和拍频信号
subplot(2,1,1);
plot(t, signal1, 'b', 'LineWidth', 1);
hold on;
plot(t, signal2, 'r', 'LineWidth', 1);
plot(t, combined_signal, 'k', 'LineWidth', 1.5);
title('原始信号和拍频信号');
xlabel('时间 (s)'); ylabel('幅度');
legend(['f1=', num2str(f1/1e3), 'kHz'], ['f2=', num2str(f2/1e3), 'kHz'], '拍频信号');
grid on;
if exist('step_idx', 'var')
    xline(t(step_idx), '--', '0.5ms', 'LabelVerticalAlignment', 'middle');
end

% 子图2：理论拍频相位
subplot(2,1,2);
plot(t, beat_phase, 'm', 'LineWidth', 1.5);
title('理论拍频相位变化');
xlabel('时间 (s)'); ylabel('相位 (rad)');
grid on;
if exist('step_idx', 'var')
    xline(t(step_idx), '--', '0.5ms', 'LabelVerticalAlignment', 'middle');
end

%% 输出拍频频率和相位斜率
beat_freq = abs(f1 - f2);
disp(['拍频频率: ', num2str(beat_freq), ' Hz']);
disp(['理论相位斜率: ', num2str(pi*(f1-f2)), ' rad/s']);
