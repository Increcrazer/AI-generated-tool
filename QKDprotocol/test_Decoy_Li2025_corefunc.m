% Li2025 成码率仿真示例；运行此脚本即可绘图。
% 输入命名、距离循环和绘图风格仿照 Lim2014/test_Decoy_Lim2014_corefunc。
addpath(fileparts(mfilename('fullpath')));

epsilon_cor = 1e-15;
epsilon_sec = 1e-15;
qX = 0.5;
N = 1e12;                     % Alice 总发送脉冲数（不是接收计数）
f = 100e6;                    % Hz
gate_width = 1e-9;            % s
width_3dB = 70e-12;           % s
k = [0.6 0.1 0.0002];         % signal, decoy, vacuum
pk = [0.5 0.25 0.25];
L = 0:2:200;                  % km
alpha = 0.2;                  % dB/km；自由空间也可直接令 alpha*L=总损耗dB
eta_Bob_detect = 0.1;
e_mis_X = 0.005;
e_mis_Z = 0.005;
f_EC = 1.16;
p_ap = 0;                     % 超导探测器一般设置0；这里不会像原Lim脚本那样忽略此输入
dc_count = [100 100 100 100];  % /s per detector
deadtime = 0;
options = struct('intensity_delta',0);

R_bitperpulse = zeros(size(L));
R_bitpersecond = zeros(size(L));
e_obs = zeros(size(L));
phi_X = zeros(size(L));
phi_Z = zeros(size(L));
nX = zeros(size(L));
details = cell(size(L));
for i = 1:numel(L)
    [R_bitperpulse(i),R_bitpersecond(i),e_obs(i),phi_X(i),nX(i),phi_Z(i),details{i}] = ...
        Decoy_Li2025_corefunc(epsilon_cor,epsilon_sec, ...
        qX,N,f,gate_width,width_3dB,k,pk,L(i),alpha,eta_Bob_detect, ...
        e_mis_X,e_mis_Z,f_EC,p_ap,dc_count,deadtime,options);
end

fprintf('First distance: %.0f km, key=%g bit, rate=%.6g bit/pulse (%.6g bit/s)\n', ...
    L(1),details{1}.l,R_bitperpulse(1),R_bitpersecond(1));
fprintf('Security: epsilon_sec=%.3g, xi=%.3g, total epsilon=%.3g\n', ...
    details{1}.security.epsilon_sec_used,details{1}.security.xi, ...
    details{1}.security.epsilon_total);

figure('Color','white','Position',[100 100 800 600]);
semilogy(L,R_bitperpulse,'LineWidth',2.5,'Color',[0 0.4470 0.7410]);
grid on;
xlabel('Distance (km)');
ylabel('Secret key rate (bit/pulse)');
title('Li 2025 finite-key analysis: signal states in X and Z');
set(gca,'FontSize',12,'LineWidth',1.5,'YMinorTick','on','XMinorTick','on');

% 实测使用方式（须改成自己当前块的计数）：
% measured = struct('n_X',[Nsx,Ndx,N0x], 'n_Z',[Nsz,Ndz,N0z], ...
%                   'm_X',[Msx,Mdx,M0x], 'm_Z',[Msz,Mdz,M0z]);
% options.counts = measured;
% options.lambda_EC = dropLen;  % 两基信号态总纠错泄露；不包括已单独减去的安全项
% [R, R_Hz, eX, phiX, nX, phiZ, result] = Decoy_Li2025_corefunc( ...
%     epsilon_cor,epsilon_sec,qX,N,f,gate_width,width_3dB,k,pk, ...
%     0,1,eta_Bob_detect,e_mis_X,e_mis_Z,f_EC,p_ap,dc_count,deadtime,options);
