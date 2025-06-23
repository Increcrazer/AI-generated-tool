% 迭代收敛到 pi 的 MATLAB 实现
rng(45);  % 设置随机种子以便复现
N_iter = 10;  % 迭代次数

% 初始化 (0~2pi之间的随机数)
phi1 = 2*pi*rand();
phi2 = 2*pi*rand();
phi3 = 0;  % 将在第一次迭代中更新

% 存储误差结果
err_phi1 = zeros(N_iter, 1);
err_phi2 = zeros(N_iter, 1);

fprintf('初始值: phi1_0 = %.4f, phi2_0 = %.4f (与π的误差: %.4f, %.4f)\n', ...
        phi1, phi2, abs(phi1-pi), abs(phi2-pi));

for n = 1:N_iter
    % 第一步: 更新 phi3
    base = (pi - (phi1 - phi2)) / 2;
    phi3 = mod(base, 2*pi);
%     if abs(mod(phi1-2*pi,2*pi))<0.5|| abs(mod(phi1-2*pi,2*pi))>6.28-0.5
%         phi3 = mod(base+pi, 2*pi);
%     end
    old_phi1 = phi1;
    % 第二步: 更新 phi1
    delta1 = phi3 - phi2/2;
    numerator1 = -2 * cos(phi2/2) * sin(delta1);
    denominator1 = 1 + 2 * cos(phi2/2) * cos(delta1);
    phi1 = atan2(numerator1, denominator1) + pi;
    
    % 第三步: 更新 phi2
    delta2 = phi3 + old_phi1/2;
    numerator2 = 2 * cos(old_phi1/2) * sin(delta2);
    denominator2 = 1 + 2 * cos(old_phi1/2) * cos(delta2);
    phi2 = atan2(numerator2, denominator2) + pi;
    
    % 计算与π的绝对误差
    err_phi1(n) = abs(phi1 - pi);
    err_phi2(n) = abs(phi2 - pi);
    
    fprintf('迭代 %d: phi1 = %.6f (误差: %.6f), phi2 = %.6f (误差: %.6f)\n', ...
            n, phi1, err_phi1(n), phi2, err_phi2(n));
end

% 绘图展示收敛过程
figure;
semilogy(1:N_iter, err_phi1, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '\phi_1 误差');
hold on;
semilogy(1:N_iter, err_phi2, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '\phi_2 误差');

grid on;
xlabel('迭代次数 n');
ylabel('与π的绝对误差 (log scale)');
title('\phi_1^n 和 \phi_2^n 收敛到 \pi 的过程');
legend('Location', 'best');
set(gca, 'FontSize', 12, 'YScale', 'log');
xlim([1 N_iter]);

% 添加参考线展示指数衰减
hold on;
ref_line = 2.^(-2.^(1:N_iter));  % 指数衰减参考线
semilogy(1:N_iter, ref_line, 'k--', 'DisplayName', '超指数衰减参考线');
