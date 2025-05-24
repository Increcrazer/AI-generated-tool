% 发单态周期，根据调平衡的两路超导通道计数反算Bloch向量夹角，并绘制调制电压-夹角曲线

% 读取Excel数据
data = readtable('1p25G period.xlsx', 'Range', 'A1:G15', 'VariableNamingRule', 'preserve');

% 提取数据列
V_rf = data.("射频电压(V)")(1:end-1);   % 射频电压列
CH3 = data.CH3(1:end-1);               % CH3数据列
CH4 = data.CH4(1:end-1);               % CH4数据列
V_dc = data.("直流电压(V)")(1:end-1);   % 直流电压列
angle_data = data.("角度")(1:end-1);        % 角度列

% 计算P = CH3/CH4
P = CH3 ./ CH4;

% 计算theta值
theta_rad = 2 * atan(sqrt(1./P));      % 计算弧度值
theta_deg = rad2deg(theta_rad);        % 转换为角度
theta_deg = theta_deg - theta_deg(1);  % 以0V为基准

% 创建图形
figure;
hold on;

% 绘制射频电压对应的角度（蓝色曲线）
plot(V_rf, theta_deg, 'bx-', 'LineWidth', 1.5);

% 找出有直流电压数据的点
dc_idx = find(~isnan(V_dc));
V_dc_valid = V_dc(dc_idx);
angle_dc_valid = angle_data(dc_idx);

% 绘制直流电压对应的角度（红色点和连线）
plot(V_dc_valid, angle_dc_valid, 'ro-', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

% 标记数据点
text(V_rf, theta_deg, num2str(theta_deg, '%.1f°'), 'VerticalAlignment','bottom');
text(V_dc_valid, angle_dc_valid, num2str(angle_dc_valid, '%.1f°'), 'Color','red', 'VerticalAlignment','top');

% 图形标注
xlabel('电压 (V)');
ylabel('角度 (°)');
title('电压-角度关系曲线');
legend('射频电压角度', '直流电压角度', 'Location', 'best');
grid on;

% 设置坐标轴范围
xlim([min(V_rf)-0.1, max(V_rf)+0.1]);
ylim([min([theta_deg; angle_dc_valid])-10, max([theta_deg; angle_dc_valid])+10]);

hold off;
