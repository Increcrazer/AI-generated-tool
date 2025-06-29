plot_phi_surfaces();

seed = 42;
rng(seed);
phi1_init = 2*pi*rand();
phi2_init = 2*pi*rand();
phi3_init = 2*pi*rand();
plot_one_iterations(phi1_init, phi2_init, phi3_init, 6);

plot_nine_iterations();

%% 对于0~2pi之间的相位，绘制单次迭代后的输出
function plot_phi_surfaces()
    % 定义输入范围 [0, 2π]
    phi1_prev_range = linspace(0, 2*pi, 100);
    phi2_prev_range = linspace(0, 2*pi, 100);
    [phi1_grid, phi2_grid] = meshgrid(phi1_prev_range, phi2_prev_range);
    
    % 预分配结果矩阵
    phi1_n_results = zeros(size(phi1_grid));
    phi2_n_results = zeros(size(phi2_grid));
    
    % 计算所有点的结果
    for i = 1:numel(phi1_grid)
        [phi1_n_results(i), phi2_n_results(i),~,~] = calculate_phis(phi1_grid(i), phi2_grid(i),0,0);
    end
    
    % 创建图形窗口
    figure('Position', [100, 100, 1200, 500]);
    
    % 绘制 phi1^n 曲面
    subplot(1, 2, 1);
    surf(phi1_grid, phi2_grid, phi1_n_results, 'EdgeColor', 'none');
    xlabel('\phi_1^{n-1}');
    ylabel('\phi_2^{n-1}');
    zlabel('\phi_1^n');
    title('\phi_1^n 曲面');
    colormap('jet');
    colorbar;
    view(3);
    axis tight;
    
    % 绘制 phi2^n 曲面
    subplot(1, 2, 2);
    surf(phi1_grid, phi2_grid, phi2_n_results, 'EdgeColor', 'none');
    xlabel('\phi_1^{n-1}');
    ylabel('\phi_2^{n-1}');
    zlabel('\phi_2^n');
    title('\phi_2^n 曲面');
    colormap('jet');
    colorbar;
    view(3);
    axis tight;
end

%% 根据九种初始相位，计算每一次迭代后的相位，并绘图
function plot_nine_iterations()
    % 创建3×3的子图布局
    figure('Position', [100, 100, 1400, 1400], 'Color', 'white');
    
    % 定义9个固定种子（确保可复现性）
    seeds = [42, 123, 7, 1984, 2023, 314, 2718, 1618, 100];
    
    for i = 1:9
        % 设置当前子图
        subplot(3, 3, i);
        
        % 固定随机种子并生成初始值
        rng(seeds(i));
        phi1_init = 2*pi*rand();
        phi2_init = 2*pi*rand();
        phi3_init = 2*pi*rand();
        
        % 模拟迭代数据
        num_iterations = 6;
        [phi1_list, phi2_list, phi3_list] = iterate_phis(phi1_init, phi2_init, phi3_init, num_iterations);
        phi1_err_list = abs(phi1_list - pi);
        phi2_err_list = abs(phi2_list - pi);
        phi3_err_list = min(abs(phi3_list - 1/2*pi), abs(phi3_list - 3/2*pi));
        
        % ========== 左侧坐标轴：原始相位角 ==========
        yyaxis left;
        
        % 绘制三条相位曲线
        p1 = plot(0:num_iterations, phi1_list, 'Color', [0, 0.45, 0.74], ...
                  'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', ...
                  'MarkerSize', 4, 'MarkerFaceColor', [0, 0.45, 0.74], ...
                  'DisplayName', '\phi_1^n');
        hold on;
        p2 = plot(0:num_iterations, phi2_list, 'Color', [0.85, 0.33, 0.1], ...
                  'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'square', ...
                  'MarkerSize', 4, 'MarkerFaceColor', [0.85, 0.33, 0.1], ...
                  'DisplayName', '\phi_2^n');

        p3 = plot(0:num_iterations, phi3_list, 'Color', [0.47, 0.67, 0.19], ...
                  'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'diamond', ...
                  'MarkerSize', 4, 'MarkerFaceColor', [0.47, 0.67, 0.19], ...
                  'DisplayName', '\phi_3^n');
        
        % 左侧坐标轴设置
        ax = gca;
        ax.YColor = 'k';
        ax.FontSize = 9;
        xlim([0 num_iterations]);
        ylim([0 2*pi]);
        xticks(0:2:num_iterations);
        yticks(0:pi/2:2*pi);
        yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
        grid on;
        
        % 添加标题（显示种子）
         title(sprintf('\\phi_1^0=%.2f, \\phi_2^0=%.2f',...
              phi1_init, phi2_init), 'FontSize', 10);
        
        % ========== 右侧坐标轴：误差对数坐标 ==========
        yyaxis right;
        
        % 绘制误差曲线（虚线样式）
        e1 = semilogy(0:num_iterations, phi1_err_list, ':', 'Color', [0, 0.45, 0.74], ...
                     'LineWidth', 1.5, 'DisplayName', '|\phi_1^n-\pi|');
        hold on;
        e2 = semilogy(0:num_iterations, phi2_err_list, ':', 'Color', [0.85, 0.33, 0.1], ...
                     'LineWidth', 1.5, 'DisplayName', '|\phi_2^n-\pi|');
        e3 = semilogy(0:num_iterations, phi3_err_list, ':', 'Color', [0.47, 0.67, 0.19], ...
                     'LineWidth', 1.5, 'DisplayName', 'min|\phi_3^n-\pi/2|');
        
        % 添加超指数衰减参考线
        n = 0:num_iterations;
        ref_line = pi*exp(-exp(0.8*n));  % 可调整参数
        r = semilogy(n, ref_line, 'k--', 'LineWidth', 1);
        
        % 右侧坐标轴设置
        ax.YColor = 'k';
        ax.YScale = 'log';
        ylim([1e-6, pi]);  % 根据实际误差调整
        
        % 添加坐标轴标签（仅边缘子图显示）
        if ismember(i, [1,4,7])
            yyaxis left;
            ylabel('Angle (rad)', 'FontSize', 9);
            yyaxis right;
            ylabel('Error (log scale)', 'FontSize', 9);
        end
        if ismember(i, [7,8,9])
            xlabel('Iteration number', 'FontSize', 9);
        end
    end
    
    % 添加共享图例（放置在图形底部中央）
    leg = legend([p1, p2, p3, e1, e2, e3, r], ...
                {'\phi_1^n', '\phi_2^n', '\phi_3^n', ...
                 '|\phi_1^n-\pi|', '|\phi_2^n-\pi|', '|\phi_3^n-1(3)\pi/2|', ...
                 'Super-exp'}, ...
                'Orientation', 'horizontal', ...
                'NumColumns', 3, ...
                'Position', [0.65 0.02 0.15 0.03]);
    leg.FontSize = 9;
    leg.Box = 'off';
end
%% 根据一种初始相位，计算每一次迭代后的相位，并绘图
function plot_one_iterations(phi1_init, phi2_init, phi3_init, num_iterations)
    [phi1_list, phi2_list, phi3_list] = iterate_phis(phi1_init, phi2_init, phi3_init, num_iterations);
    phi1_err_list = abs(phi1_list - pi);
    phi2_err_list = abs(phi2_list - pi);
    phi3_err_list = min(abs(phi3_list - 1/2*pi), abs(phi3_list - 3/2*pi));
    
    % 创建图形窗口
    figure('Position', [100, 100, 1000, 600], 'Color', 'white');
    
    % ========== 左侧坐标轴：原始相位角 ==========
    yyaxis left;
    
    % 绘制三条相位曲线
    p1 = plot(0:num_iterations, phi1_list, 'Color', [0, 0.45, 0.74], ...
              'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', [0, 0.45, 0.74], ...
              'DisplayName', '\phi_1^n');
    hold on;
    
    p2 = plot(0:num_iterations, phi2_list, 'Color', [0.85, 0.33, 0.1], ...
              'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'square', ...
              'MarkerSize', 8, 'MarkerFaceColor', [0.85, 0.33, 0.1], ...
              'DisplayName', '\phi_2^n');
    
    p3 = plot(0:num_iterations, phi3_list, 'Color', [0.47, 0.67, 0.19], ...
              'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'diamond', ...
              'MarkerSize', 8, 'MarkerFaceColor', [0.47, 0.67, 0.19], ...
              'DisplayName', '\phi_3^n');
    
    % 左侧坐标轴设置
    ax = gca;
    ax.YColor = 'k';
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    ax.LineWidth = 1.2;
    ax.XLabel.String = "Iteration number";
    ax.XLabel.FontSize = 14;
    ax.YLabel.String = "Phase angle (rad)";
    ax.YLabel.FontSize = 14;
    ax.Title.String = sprintf('Iteration Process (\\phi_1^0=%.2f, \\phi_2^0=%.2f)', ...
        phi1_init, phi2_init);
    ax.Title.FontSize = 16;
    ax.Title.FontWeight = 'bold';
    
    xlim([0, num_iterations]);
    ylim([0, 2*pi]);
    xticks(0:2:num_iterations);
    yticks(0:0.5*pi:2*pi);
    yticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    grid on;
    box on;
    
    % ========== 右侧坐标轴：误差对数坐标 ==========
    yyaxis right;
    
    % 绘制误差曲线（虚线样式）
    e1 = semilogy(0:num_iterations, phi1_err_list, ':', 'Color', [0, 0.45, 0.74], ...
                 'LineWidth', 2, 'DisplayName', '|\phi_1^n-\pi|');
    hold on;
    e2 = semilogy(0:num_iterations, phi2_err_list, ':', 'Color', [0.85, 0.33, 0.1], ...
                 'LineWidth', 2, 'DisplayName', '|\phi_2^n-\pi|');
    e3 = semilogy(0:num_iterations, phi3_err_list, ':', 'Color', [0.47, 0.67, 0.19], ...
                 'LineWidth', 2, 'DisplayName', 'min|\phi_3^n-\pi/2|');
    
    % 添加超指数衰减参考线 (a*exp(-b*exp(c*n)))
    n = 0:num_iterations;
    ref_line = 2*exp(-exp(0.5*n));  % 可调整参数
    r = semilogy(n, ref_line, 'k--', 'LineWidth', 1.5, ...
                'DisplayName', 'Super-exp decay');
    
    % 右侧坐标轴设置
    ax.YColor = 'k';
    ax.YLabel.String = 'Error (log scale)';
    ax.YLabel.FontSize = 14;
    ax.YScale = 'log';
    ylim([1e-6, pi]);  % 根据实际误差范围调整
    
    % ========== 组合图例 ==========
    leg = legend([p1, p2, p3, e1, e2, e3, r], ...
                'Location', 'northeast', 'NumColumns', 1);
    leg.FontSize = 10;
    leg.Box = 'off';
    
end

%% 根据迭代次数和初始相位，计算每一次迭代后的相位
function [phi1_list, phi2_list, phi3_list] = iterate_phis(phi1_init, phi2_init, phi3_init, num_iterations)
    % 初始化列表（包括初始值）
    phi1_list = zeros(1, num_iterations + 1);
    phi2_list = zeros(1, num_iterations + 1);
    phi3_list = zeros(1, num_iterations + 1);
    I_list = zeros(1, num_iterations + 1);
    
    % 存储初始值
    phi1_list(1) = mod(phi1_init, 2*pi);
    phi2_list(1) = mod(phi2_init, 2*pi);
    phi3_list(1) = mod(phi3_init, 2*pi);
    
    % 迭代计算
    for n = 1:num_iterations
        [phi1_list(n+1), phi2_list(n+1), phi3_list(n+1), I_list(n+1)] = calculate_phis(phi1_list(n), phi2_list(n), phi3_list(n),0);
        if I_list(n+1)>I_list(n) && n>1
            phi3_change = mod(phi3_list(n) + pi,2*pi);
            [phi1_list(n+1), phi2_list(n+1), phi3_list(n+1), I_list(n+1)] = calculate_phis(phi1_list(n), phi2_list(n), phi3_change,1);
        end
    end
end

%% 迭代算法表达式
function [phi1_n, phi2_n, phi3_n, I_n] = calculate_phis(phi1_n_minus_1, phi2_n_minus_1, phi3_n_minus_1, flag)
    % 第一步：计算 phi3_n以及该步的光强
    k = 0; % 可以调整k值
    
    if flag == 0
        phi3_n_star = (k + 0.5)*pi - 0.5*(phi1_n_minus_1 - phi2_n_minus_1);
        phi3_n = mod(phi3_n_star, 2*pi);
    elseif flag ==1
        phi3_n = mod(phi3_n_minus_1, 2*pi);
    end
    I_n = 1/8*(cos(phi1_n_minus_1)+1) + 1/8*(cos(phi2_n_minus_1)+1) + ...
        1/2*cos(phi1_n_minus_1/2)*cos(phi2_n_minus_1/2) ...
        *cos((phi1_n_minus_1 - phi2_n_minus_1)/2 + phi3_n);
    
    % 第二步：计算 phi1_n
    delta = phi3_n - phi2_n_minus_1/2;
    numerator = 2 * cos(phi2_n_minus_1/2) * sin(delta);
    denominator = 1 + 2 * cos(phi2_n_minus_1/2) * cos(delta);
    phi1_n = -atan2(numerator, denominator) + pi + 2*k*pi;
    phi1_n = mod(phi1_n, 2*pi);
    
    % 第三步：计算 phi2_n
    delta = phi3_n + phi1_n/2;
    numerator = 2 * cos(phi1_n/2) * sin(delta);
    denominator = 1 + 2 * cos(phi1_n/2) * cos(delta);
    phi2_n = atan2(numerator, denominator) + pi + 2*k*pi;
    phi2_n = mod(phi2_n, 2*pi);
end
