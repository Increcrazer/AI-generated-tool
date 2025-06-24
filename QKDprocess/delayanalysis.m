%% 读取表格文件
data = readtable('1250M_PM.csv');  
time_full = (data{1430:1485, 1} - data{1430, 1}) * 1e12;  % 转换为ps单位
voltage_full = data{1430:1485, 2} * 31.6228;              % 第二列数据（调整幅度）
%% 绘制原始电压波形（figure 1）并标记200-400ps区间
figure(1);
clf;  % 清空图形
set(gcf, 'Position', [100 100 800 500]);  % 设置图形大小

% 绘制主电压波形
plot(time_full, voltage_full, 'b-', 'LineWidth', 1.5);
hold on;

% 标记200-400ps区间
x_range = [200, 400];  % 要标记的区间
y_range = [min(voltage_full), max(voltage_full)];

% 绘制区间左右边界线
plot([x_range(1) x_range(1)], y_range, 'r--', 'LineWidth', 1.2);
plot([x_range(2) x_range(2)], y_range, 'r--', 'LineWidth', 1.2);

% 添加区间填充色（半透明）
fill([x_range(1) x_range(2) x_range(2) x_range(1)], ...
     [y_range(1) y_range(1) y_range(2) y_range(2)], ...
     [1 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% % 添加文本标注
% text(mean(x_range), max(voltage_full)*0.95, ...
%     '200-400ps', ...
%     'HorizontalAlignment', 'center', ...
%     'FontSize', 11, 'Color', 'b', 'FontWeight', 'bold');

hold off;

% 图形美化
xlabel('Time (ps)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
title('Typical Voltage Waveform', 'FontSize', 16, 'FontWeight','bold');
grid on;
xlim([min(time_full), max(time_full)]);  % 横坐标顶格
% legend('Typical Voltage Waveform', 'Location', 'best');

%% 创建高密度时间网格（插值）
time_zoom = (data{1445:1465, 1} - data{1440, 1}) * 1e12;  % 原始时间点(ps)
voltage_zoom = data{1445:1465, 2} * 31.6228;              % 原始电压数据(mV)
num_points = 500;  % 增加至500个点
time_dense = linspace(min(time_zoom), max(time_zoom), num_points)';
voltage_dense = interp1(time_zoom, voltage_zoom, time_dense, 'spline');
%% 绘制图形（高分辨率版本）
figure(2);
clf;
set(gcf, 'Position', [100 100 1000 700], 'Color', 'w');

% ================= 左侧坐标：电压波形 =================
yyaxis left;
h_voltage = plot(time_dense, voltage_dense, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', 'Voltage');
ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
voltage_range = max(voltage_dense) - min(voltage_dense);
ylim([min(voltage_dense)-0.1*voltage_range, max(voltage_dense)+0.15*voltage_range]);

% ================= 右侧坐标：高斯脉冲 =================
yyaxis right;
hold on;

% 高斯脉冲参数
sigma = 60/(2*sqrt(2*log(2)));  % σ ≈ 25.5ps (FWHM=60ps)
pulse_positions = [180, 190, 200];
colors = {[0.8 0 0], [0 0.6 0], [0.9 0.4 0]};  % 深红、深绿、橙色
pulse_handles = gobjects(1,3);

% 超高分辨率高斯脉冲计算
for i = 1:3
    t0 = pulse_positions(i);
    
    % 在密集网格上计算高斯脉冲
    pulse = exp(-(time_dense-t0).^2 / (2*sigma^2));
    pulse_handles(i) = plot(time_dense, pulse, '-', 'LineWidth', 2.2, ...
        'Color', colors{i}, 'DisplayName', sprintf('at %dps(σ=%.1fps)', t0, sigma));
    
    % 3σ标记优化
    sigma3 = [t0-3*sigma, t0+3*sigma];
    plot(sigma3, [0.92 0.92], 'k:', 'LineWidth', 1.3, 'Color', [0.3 0.3 0.3]);
    plot([sigma3(1) sigma3(1)], [0 1], ':', 'LineWidth', 0.9, 'Color', [0.4 0.4 0.4]);
    plot([sigma3(2) sigma3(2)], [0 1], ':', 'LineWidth', 0.9, 'Color', [0.4 0.4 0.4]);
    
    % 改进的文本标注
    text(sigma3(1), 0.96, sprintf('%.0f', sigma3(1)), ...
        'HorizontalAlignment', 'right', 'FontSize', 11, 'Color', [0.2 0.2 0.2]);
    text(sigma3(2), 0.96, sprintf('%.0f', sigma3(2)), ...
        'HorizontalAlignment', 'left', 'FontSize', 11, 'Color', [0.2 0.2 0.2]);
end

hold off;
ylabel('Normalized Amplitude', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 1.05]);

% ================= 专业级图形美化 =================
xlabel('Time (ps)', 'FontSize', 14, 'FontWeight', 'bold');
title('Voltage & Optical pulse', ...
    'FontSize', 16, 'FontWeight', 'bold');

% 智能图例和网格
legend([h_voltage, pulse_handles], 'Location', 'southeast', 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3, 'LineWidth', 1.2);

% 优化坐标轴范围
x_span = max(time_dense) - min(time_dense);
xlim([min(time_dense)-0.05*x_span, max(time_dense)+0.05*x_span]);

% % 添加专业标注
% annotation('textbox', [0.16,0.16,0.25,0.1], 'String', ...
%     sprintf('Gaussian Parameters:\nFWHM = 60ps\nσ = %.1f ps', sigma), ...
%     'FontSize', 11, 'FitBoxToText', 'on', 'EdgeColor', 'none', ...
%     'BackgroundColor', [1 1 1 0.7]);

% 启用抗锯齿
set(gcf, 'GraphicsSmoothing', 'on');

%% 计算保真度、误码率
voltage_offset = -6.5; % 电压偏移量
voltage_dense = voltage_dense - voltage_offset; % 电压校准
V_pi = 3.0;  % 半波电压(V)
k = pi/V_pi;  % 相位调制系数

% 设置检测基矢参数（根据实际需求调整）
V0 = -2.2 - voltage_offset;  % 基准电压
phi0 = k * V0;  % 基准相位

% 初始化存储结果
pulse_ber = zeros(1,3);  
pulse_fidelity = zeros(1,3);  

for i = 1:3
    t0 = pulse_positions(i);
    
    % 1. 生成高斯脉冲并限制在±3σ范围内
    pulse = exp(-(time_dense-t0).^2 / (2*sigma^2));
    valid_idx = (time_dense >= t0-3*sigma) & (time_dense <= t0+3*sigma);
    t_valid = time_dense(valid_idx);
    pulse_valid = pulse(valid_idx);
    
    % 2. 面积归一化（∫f(t)dt=1）
    pulse_valid = pulse_valid / trapz(t_valid, pulse_valid); 
    
    % 3. 计算相位调制 ϕ(t) = k*(V(t)-V0)
    V_valid = voltage_dense(valid_idx);
    phi = k*(V_valid - V0);
    
    % 4. 计算严格归一化因子N（通过双积分）
    N = 0;
    for m = 1:length(t_valid)
        % 计算dt1（MATLAB需用if-else代替三元运算符）
        if m > 1
            dt1 = t_valid(m) - t_valid(m-1);
        else
            dt1 = t_valid(2) - t_valid(1);
        end
        
        for n = 1:length(t_valid)
            % 计算dt2
            if n > 1
                dt2 = t_valid(n) - t_valid(n-1);
            else
                dt2 = t_valid(2) - t_valid(1);
            end
            
            integrand = sqrt(pulse_valid(m)*pulse_valid(n)) * ...
                       (1 + exp(1i*(phi(n)-phi(m))))/2;
            N = N + integrand * dt1 * dt2;
        end
    end
    
    % 5. 向量化优化版本（替代双重循环，更高效）
    % [T1,T2] = meshgrid(t_valid);
    % [P1,P2] = meshgrid(pulse_valid);
    % [Phi1,Phi2] = meshgrid(phi);
    % integrand_matrix = sqrt(P1.*P2).*(1+exp(1i*(Phi1-Phi2)))/2;
    % N = trapz(t_valid, trapz(t_valid, integrand_matrix, 2));
    
    % 6. 计算未归一化的量子态重叠
    integrand_psi = sqrt(pulse_valid) .* (1 + exp(1i*phi)) /2;
    overlap_unnormalized = trapz(t_valid, integrand_psi);
    
    % 7. 最终结果
    pulse_fidelity(i) = abs(overlap_unnormalized)^2 / N;
    pulse_ber(i) = 1 - pulse_fidelity(i);
    
    fprintf('脉冲@%dps: N=%.6f, 保真度=%.6f\n', t0, N, pulse_fidelity(i));
end

%% 结果展示
fprintf('\n=== 量子误码率分析结果 ===\n');
fprintf('系统参数验证:\n');
fprintf('• 电压偏移量 = %.2f V\n', voltage_offset);
fprintf('• 半波电压 Vπ = %.2f V\n', V_pi);
fprintf('• 基准电压 V0 = %.2f V (φ0=%.2fπ)\n', V0, phi0/pi);
fprintf('• 电压范围 = [%.2f, %.2f] V\n', min(voltage_dense), max(voltage_dense));
fprintf('----------------------------------------\n');

for i = 1:3
    fprintf('脉冲@%dps:\n', pulse_positions(i));
    fprintf('• 保真度 = %.6f (合理范围验证: %s)\n', ...
            pulse_fidelity(i), iff(pulse_fidelity(i) >= 0 && pulse_fidelity(i) <= 1, '√', '×'));
    fprintf('• 误码率 = %.2e\n', pulse_ber(i));
end
fprintf('========================================\n');


%% 辅助函数
function s = iff(condition, true_str, false_str)
    if condition
        s = true_str;
    else
        s = false_str;
    end
end

