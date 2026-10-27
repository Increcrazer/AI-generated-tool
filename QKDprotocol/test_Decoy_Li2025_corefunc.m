% Li2025成码率仿真。将MATLAB当前文件夹切换到本目录后直接运行。
% 主体布局与test_Decoy_Lim2014_corefunc一致；不自动设置目录或搜索路径。
% 单位：距离km，损耗dB/km，频率Hz，时间s，暗计数/s。
% Li2025只使用X/Z基信号态成钥；以下数值为演示参数。

%% 输入参数
epsilon_cor = 1e-15;
epsilon_sec = 1e-15;
qX = 0.5;
N = 1e12;                       % Alice总发送脉冲数，不是筛后块长
f = 100e6;
gate_width = 1e-9;               % 软件筛选窗宽度
width_3dB = 70e-12;
k = [0.6 0.1 0.0002];            % 信号、诱骗、真空
pk = [0.5 0.25 0.25];
L = 0:2:200;
alpha = 0.2;
eta_Bob_detect = 0.1;
e_mis_X = 0.005;
e_mis_Z = 0.005;
f_EC = 1.16;
p_ap = 0;                       % 当前忽略SNSPD后脉冲
dc_count = [100 100 100 100];     % 各通道暗计数率，演示值
deadtime = 0;                   % 每基共用的有效死时间，演示值

%% 逐距离计算
R_bitperpulse = zeros(1,length(L));
R_bitpersecond = zeros(1,length(L));
e_obs = zeros(1,length(L));
phi_X = zeros(1,length(L));
phi_Z = zeros(1,length(L));
nX = zeros(1,length(L));

for i = 1:length(L)
    [R_bitperpulse(i), R_bitpersecond(i), e_obs(i), phi_X(i), nX(i), phi_Z(i)] = ...
        Decoy_Li2025_corefunc ...
            (epsilon_cor, epsilon_sec, ...
             qX, N, f, gate_width, width_3dB, ...
             k, pk, ...
             L(i), alpha, eta_Bob_detect, ...
             e_mis_X, e_mis_Z, f_EC, ...
             p_ap, dc_count, deadtime);
end

%% 绘图
figure('Color','white','Position',[100 100 800 600]);
semilogy(L,R_bitperpulse,'LineWidth',2.5,'Color',[0 0.4470 0.7410]);
grid on;
xlabel('Distance (km)','FontSize',14,'FontWeight','bold');
ylabel('Secret key rate (bit/pulse)','FontSize',14,'FontWeight','bold');
title('Li2025 key rate vs distance','FontSize',16,'FontWeight','bold');

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.YScale = 'log';
ax.YMinorTick = 'on';
ax.XMinorTick = 'on';

% 光强波动、自定义安全参数和实测计数仍可通过第19个输入options指定。
% 第7个输出detail可查看中间统计量；具体调用示例见README.md。
