%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023/05 @copyright WYA
% IM和IQ调制器的完全体建模，包括PS、PDL、BS不完美性
% 仿真PIC中IQ和IM在不完美BS情况下的消光比
% IM和IQ调制器电压各考虑两种方案
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Part
%% Test validity of the function
Ein = 1;
R = 0.5;

k = -0.8;
b = 0;
CDM_L = 0.6;    % 6 mm
PS_L = 300e-4;  % 300 um
prop_loss = 0;  % 3 dB/cm

%% define IM and IQ working points
params = struct();

% % IM
% params.CDM_fai1_IM_S = pi;
% params.CDM_fai2_IM_S = 0;
% params.PS_fai1_IM_S = 0;
% params.PS_fai2_IM_S = pi; 
% params.CDM_fai1_IM_D = -acos(-1/2)+pi;
% params.CDM_fai2_IM_D = 0;
% params.PS_fai1_IM_D = 0;
% params.PS_fai2_IM_D = pi;
% params.CDM_fai1_IM_V = 0;
% params.CDM_fai2_IM_V = 0;
% params.PS_fai1_IM_V = 0;
% params.PS_fai2_IM_V = pi;

params.CDM_fai1_IM_S = 0;
params.CDM_fai2_IM_S = acos(-1/2);
params.PS_fai1_IM_S = acos(-1/2);
params.PS_fai2_IM_S = 0;
params.CDM_fai1_IM_D = 0;
params.CDM_fai2_IM_D = 0;
params.PS_fai1_IM_D = acos(-1/2);
params.PS_fai2_IM_D = 0;
params.CDM_fai1_IM_V = pi-acos(-1/2);
params.CDM_fai2_IM_V = 0;
params.PS_fai1_IM_V = acos(-1/2);
params.PS_fai2_IM_V = 0;

% IQ 
% params.CDM_fai2_IQ_S = pi;
% params.CDM_fai4_IQ_S = 0;
% params.PS_fai11_IQ_S = pi;
% params.PS_fai22_IQ_S = 0;
% params.PS_fai2_IQ_S = pi;
% params.CDM_fai2_IQ_D = 0;
% params.CDM_fai4_IQ_D = 0;
% params.PS_fai11_IQ_D = pi;
% params.PS_fai22_IQ_D = 0;
% params.PS_fai2_IQ_D = pi;
% params.CDM_fai2_IQ_V = 0;
% params.CDM_fai4_IQ_V = pi;
% params.PS_fai11_IQ_V = pi;
% params.PS_fai22_IQ_V = 0;
% params.PS_fai2_IQ_V = pi;

params.CDM_fai2_IQ_S = 0;
params.CDM_fai4_IQ_S = 0;
params.PS_fai11_IQ_S = 0;
params.PS_fai22_IQ_S = 0;
params.PS_fai2_IQ_S = 0;
params.CDM_fai2_IQ_D = 0;
params.CDM_fai4_IQ_D = 0;
params.PS_fai11_IQ_D = pi;
params.PS_fai22_IQ_D = 0;
params.PS_fai2_IQ_D = 0;
params.CDM_fai2_IQ_V = pi;
params.CDM_fai4_IQ_V = pi;
params.PS_fai11_IQ_V = 0;
params.PS_fai22_IQ_V = 0;
params.PS_fai2_IQ_V = pi;


%% test IM
[~,IM_S] = MZI(Ein,0,R,R,...
        params.CDM_fai1_IM_S,params.CDM_fai2_IM_S,CDM_L,k,b,...
        params.PS_fai1_IM_S,params.PS_fai2_IM_S,PS_L,...
        prop_loss);

[~,IM_D] = MZI(Ein,0,R,R,...
        params.CDM_fai1_IM_D,params.CDM_fai2_IM_D,CDM_L,k,b,...
        params.PS_fai1_IM_D,params.PS_fai2_IM_D,PS_L,...
        prop_loss);

[~,IM_V] = MZI(Ein,0,R,R,...
        params.CDM_fai1_IM_V,params.CDM_fai2_IM_V,CDM_L,k,b,...
        params.PS_fai1_IM_V,params.PS_fai2_IM_V,PS_L,...
        prop_loss);
    
disp('SDV intensity proportion of IM is:');
disp([abs(IM_S)^2/Ein^2, abs(IM_D)^2/Ein^2, abs(IM_V)^2/Ein^2]);

%% test IQ
IQ_S = layout_IQ(Ein,R,R,R,R,R,R, ...
        params.CDM_fai2_IQ_S,params.CDM_fai4_IQ_S,CDM_L,k,b,  ...
        params.PS_fai11_IQ_S,params.PS_fai22_IQ_S,params.PS_fai2_IQ_S,PS_L, ...
        prop_loss);
    
IQ_D = layout_IQ(Ein,R,R,R,R,R,R, ...
        params.CDM_fai2_IQ_D,params.CDM_fai4_IQ_D,CDM_L,k,b,  ...
        params.PS_fai11_IQ_D,params.PS_fai22_IQ_D,params.PS_fai2_IQ_D,PS_L, ...
        prop_loss);
    
IQ_V = layout_IQ(Ein,R,R,R,R,R,R, ...
        params.CDM_fai2_IQ_V,params.CDM_fai4_IQ_V,CDM_L,k,b,  ...
        params.PS_fai11_IQ_V,params.PS_fai22_IQ_V,params.PS_fai2_IQ_V,PS_L, ...
        prop_loss);

disp('SDV intensity proportion of IQ is:');
disp([abs(IQ_S)^2/Ein^2, abs(IQ_D)^2/Ein^2, abs(IQ_V)^2/Ein^2]);

%% Compare ER of IM and IQ as k,b, assuming all BS of PIC are different
k = linspace(-0.8,0,40);
b = 7.7;
CDM_L = 0.6;    %6 mm
PS_L = 300e-4;  % 300 um
prop_loss = 3;  %3 dB/cm

ER_IM = zeros(1,40);
ER_IQ = zeros(1,40);
for i = 1:40
    [ER_IM(i),ER_IQ(i)] = ER(k(i),b,CDM_L,PS_L,prop_loss,params);
end

figure(1);
plot(k, ER_IM, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6);
hold on;
plot(k, ER_IQ, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--', 'Marker', 's', 'MarkerSize', 6);
hold off;

% 美化图形
grid on; % 添加网格
box on;  % 添加边框

% 坐标轴标签
xlabel('k (Phase Loss Coefficient)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Extinction Ratio (dB)', 'FontSize', 12, 'FontWeight', 'bold');

% 标题
title('Extinction Ratio Comparison: IM vs. IQ Modulators', 'FontSize', 14, 'FontWeight', 'bold');

% 图例
legend('IM Modulator', 'IQ Modulator', 'Location', 'best', 'FontSize', 10);

% 坐标轴范围调整（可选）
xlim([min(k), max(k)]); % 确保 x 轴范围覆盖所有数据
ylim([20, 40]); % 根据实际 ER 范围调整

% 设置坐标轴字体
set(gca, 'FontSize', 10, 'FontWeight', 'bold');

%% Function Part
%% calculate ER @loss(k,b) for IM and IQ, assuming all BS of PIC are different, Monte Carlo method
function [ER_mean_IM,ER_mean_IQ] = ER(k,b,CDM_L,PS_L,prop_loss,params)
    Ein = 1;
    R1 = normrnd(1/2,0.02,1,10000);
    R2 = normrnd(1/2,0.02,1,10000);
    R3 = normrnd(1/2,0.02,1,10000);
    R4 = normrnd(1/2,0.02,1,10000);
    R5 = normrnd(1/2,0.02,1,10000);
    R6 = normrnd(1/2,0.02,1,10000);

    ER_IM = zeros(1,10000);
    ER_IQ = zeros(1,10000);
    for i = 1:10000
            [~,IM_out_S] = MZI(Ein,0,R1(i),R2(i),...
                            params.CDM_fai1_IM_S,params.CDM_fai2_IM_S,CDM_L,k,b,...
                            params.PS_fai1_IM_S,params.PS_fai2_IM_S,PS_L,...
                            prop_loss);
            [~,IM_out_V] = MZI(Ein,0,R1(i),R2(i),...
                            params.CDM_fai1_IM_V,params.CDM_fai2_IM_V,CDM_L,k,b,...
                            params.PS_fai1_IM_V,params.PS_fai2_IM_V,PS_L,...
                            prop_loss);
            IQ_out_S = layout_IQ(Ein,R1(i),R2(i),R3(i),R4(i),R5(i),R6(i), ...
                            params.CDM_fai2_IQ_S,params.CDM_fai4_IQ_S,CDM_L,k,b,  ...
                            params.PS_fai11_IQ_S,params.PS_fai22_IQ_S,params.PS_fai2_IQ_S,PS_L, ...
                            prop_loss);
            IQ_out_V = layout_IQ(Ein,R1(i),R2(i),R3(i),R4(i),R5(i),R6(i), ...
                            params.CDM_fai2_IQ_V,params.CDM_fai4_IQ_V,CDM_L,k,b,  ...
                            params.PS_fai11_IQ_V,params.PS_fai22_IQ_V,params.PS_fai2_IQ_V,PS_L, ...
                            prop_loss);
            IM_out_S = abs(IM_out_S)^2;
            IM_out_V = abs(IM_out_V)^2;
            IQ_out_S = abs(IQ_out_S)^2;
            IQ_out_V = abs(IQ_out_V)^2;
            ER_IM(i) = 10*log10(IM_out_S/IM_out_V);
            ER_IQ(i) = 10*log10(IQ_out_S/IQ_out_V);
    end
    ER_mean_IM = mean(ER_IM);
    ER_mean_IQ = mean(ER_IQ);
end

%% 2to2 IQ in layout
function Eo = layout_IQ(IQ_IN,R1,R2,R3,R4,R5,R6, ...
CDM_fai2,CDM_fai4,CDM_L,k,b,  ...
PS_fai11,PS_fai22,PS_fai2,PS_L, ...
prop_loss)
    % 调用IQ函数
    [~,Eo] = IQ(IQ_IN,0,R1,R2,R3,R4,R5,R6, ...
            0,CDM_fai2,0,CDM_fai4,CDM_L,k,b, ...
            PS_fai11,0,0,PS_fai22,0,PS_fai2,PS_L, ...
            prop_loss);
end

%% 2to2 CDM-IQ 
function [Eo1,Eo2] = IQ(Ein1,Ein2,R1,R2,R3,R4,R5,R6,...
                        CDM_fai1,CDM_fai2,CDM_fai3,CDM_fai4,CDM_L,k,b,...
                        PS_fai11,PS_fai12,PS_fai21,PS_fai22,PS_fai1,PS_fai2,PS_L,...
                        prop_loss)
    [E1,E2] = BS(R1,Ein1,Ein2);
    [~,E3] =  MZI(E1,0,R2,R3,...
                CDM_fai1,CDM_fai2,CDM_L,k,b,...
                PS_fai11,PS_fai12,PS_L,...
                prop_loss);
    [~,E4] =  MZI(E2,0,R4,R5,...
                CDM_fai3,CDM_fai4,CDM_L,k,b,...
                PS_fai21,PS_fai22,PS_L,...
                prop_loss);
    E5 = E3*exp(1i*PS_fai1)*sqrt(10^(-prop_loss*PS_L/10));
    E6 = E4*exp(1i*PS_fai2)*sqrt(10^(-prop_loss*PS_L/10));
    [Eo1,Eo2] = BS(R6,E5,E6);
end

%% 2to2 CDM-MZI,including PS
% prop_loss unit is dB/cm
% CDM_L & PS_L unit is cm
function [Eo1,Eo2] = MZI(Ein1,Ein2,R1,R2,...
                        CDM_fai1,CDM_fai2,CDM_L,k,b,...
                        PS_fai1,PS_fai2,PS_L,...
                        prop_loss)
    loss_factor1 = k*CDM_fai1+b+prop_loss*CDM_L;   % k is negative
    loss_factor2 = k*CDM_fai2+b+prop_loss*CDM_L;
    loss_factor3 = prop_loss*PS_L;
    loss_factor4 = prop_loss*PS_L;
    [E1,E2] = BS(R1,Ein1,Ein2);
    E3 = E1*exp(1i*CDM_fai1)*sqrt(10^(-loss_factor1/10));
    E4 = E2*exp(1i*CDM_fai2)*sqrt(10^(-loss_factor2/10));
    E5 = E3*exp(1i*PS_fai1)*sqrt(10^(-loss_factor3/10));
    E6 = E4*exp(1i*PS_fai2)*sqrt(10^(-loss_factor4/10));
    [Eo1,Eo2] = BS(R2,E5,E6);
end

%% 2to2 BS (MMI)
function [Eo1,Eo2] = BS(R,Ein1,Ein2)
    T = 1-R;
    Eo1 = 1i*sqrt(R)*Ein1 + sqrt(T)*Ein2;
    Eo2 = sqrt(T)*Ein1 + 1i*sqrt(R)*Ein2;
end

% % 更换需要注意改MZI和IQ输出端口
% %% 2to2 BS (Y branch)
% function [Eo1,Eo2] = BS(R,Ein1,Ein2)
%     T = 1-R;
%     Eo1 = sqrt(R)*Ein1 + sqrt(T)*Ein2;
%     Eo2 = sqrt(T)*Ein1 - sqrt(R)*Ein2;
% end
