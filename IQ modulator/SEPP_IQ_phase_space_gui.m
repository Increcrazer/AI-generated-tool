function SEPP_IQ_phase_space_gui
    % 创建主窗口（支持全屏）
    fig = figure('Name', 'IQ Modulator Simulator', 'NumberTitle', 'off', ...
                'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], ...
                'Resize', 'on', 'MenuBar', 'none', 'ToolBar', 'none');
    
    % 默认参数
    params.IQ_IN = 1;
    params.R = 0.5;
    params.IQ_fai0 = 0;
    params.up_fai0 = 0;
    params.down_fai0 = 0;
    params.I_RF_phi = 0;    % I路RF相位
    params.Q_RF_phi = 0;    % Q路RF相位
    params.k = 0;
    params.b = 0;
    params.PS2_fai = 0;     % 固定为0
    params.PS3_fai = 0;     % 固定为0
    params.PS5_fai = 0;     % 输出相位调整
    
    % 创建主面板（替代uigridlayout）
    mainPanel = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
    
    % 星座图区域（顶部70%空间）
    axPanel = uipanel('Parent', mainPanel, 'Units', 'normalized', 'Position', [0.1 0.3 0.8 0.65]);
    ax = axes('Parent', axPanel, 'OuterPosition', [0 0 1 1]);
    axis(ax, 'equal');
    xlim(ax, [-1.2 1.2]);
    ylim(ax, [-1.2 1.2]);
    grid(ax, 'on');
    title(ax, 'IQ Constellation Diagram');
    xlabel(ax, 'In-phase (I)');
    ylabel(ax, 'Quadrature (Q)');
    hold(ax, 'on');
    
    % 绘制单位圆
    theta = linspace(0, 2*pi, 100);
    plot(ax, cos(theta), sin(theta), 'k--', 'LineWidth', 0.5);
    hold(ax, 'off');
    
    % 控制面板区域（底部30%空间）
    controlPanel = uipanel('Parent', mainPanel, 'Title', 'Controls', ...
                          'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.25]);
    
    % 创建三个滑块控件（使用传统UI组件）
    createTraditionalControls(controlPanel, params, ax);
    
    % 初始绘图
    updatePlot(fig, ax, params);
end

function createTraditionalControls(parent, params, ax)
    % I_RF_phi 控件
    uicontrol('Parent', parent, 'Style', 'text', ...
              'String', 'I RF Phase (0-2π):', ...
              'Units', 'normalized', 'Position', [0.05 0.8 0.25 0.15], ...
              'HorizontalAlignment', 'left');
    
    slider1 = uicontrol('Parent', parent, 'Style', 'slider', ...
                       'Min', 0, 'Max', 6.28, 'Value', params.I_RF_phi, ...
                       'SliderStep', [0.01/(6.28-0) 0.1/(6.28-0)], ...
                       'Units', 'normalized', 'Position', [0.3 0.8 0.4 0.15], ...
                       'Callback', @(src,~) sliderCallback(src, ax, 'I_RF'));
    
    valText1 = uicontrol('Parent', parent, 'Style', 'text', ...
                        'String', '0.00', ...
                        'Units', 'normalized', 'Position', [0.75 0.8 0.2 0.15], ...
                        'HorizontalAlignment', 'left');
    set(slider1, 'UserData', valText1);
    
    % Q_RF_phi 控件
    uicontrol('Parent', parent, 'Style', 'text', ...
              'String', 'Q RF Phase (0-2π):', ...
              'Units', 'normalized', 'Position', [0.05 0.5 0.25 0.15], ...
              'HorizontalAlignment', 'left');
    
    slider2 = uicontrol('Parent', parent, 'Style', 'slider', ...
                       'Min', 0, 'Max', 6.28, 'Value', params.Q_RF_phi, ...
                       'SliderStep', [0.01/(6.28-0) 0.1/(6.28-0)], ...
                       'Units', 'normalized', 'Position', [0.3 0.5 0.4 0.15], ...
                       'Callback', @(src,~) sliderCallback(src, ax, 'Q_RF'));
    
    valText2 = uicontrol('Parent', parent, 'Style', 'text', ...
                        'String', '0.00', ...
                        'Units', 'normalized', 'Position', [0.75 0.5 0.2 0.15], ...
                        'HorizontalAlignment', 'left');
    set(slider2, 'UserData', valText2);
    
    % PS5_fai 控件
    uicontrol('Parent', parent, 'Style', 'text', ...
              'String', 'Output Phase Shift (0-2π):', ...
              'Units', 'normalized', 'Position', [0.05 0.2 0.25 0.15], ...
              'HorizontalAlignment', 'left');
    
    slider3 = uicontrol('Parent', parent, 'Style', 'slider', ...
                       'Min', 0, 'Max', 6.28, 'Value', params.PS5_fai, ...
                       'SliderStep', [0.01/(6.28-0) 0.1/(6.28-0)], ...
                       'Units', 'normalized', 'Position', [0.3 0.2 0.4 0.15], ...
                       'Callback', @(src,~) sliderCallback(src, ax, 'PS5'));
    
    valText3 = uicontrol('Parent', parent, 'Style', 'text', ...
                        'String', '0.00', ...
                        'Units', 'normalized', 'Position', [0.75 0.2 0.2 0.15], ...
                        'HorizontalAlignment', 'left');
    set(slider3, 'UserData', valText3);
    
    % 信息显示区域
    infoPanel = uipanel('Parent', parent, 'Title', 'Current Values', ...
                       'Units', 'normalized', 'Position', [0.85 0.05 0.1 0.9]);
    
    uicontrol('Parent', infoPanel, 'Style', 'text', ...
             'String', 'I: 0.00\nQ: 0.00\nAmplitude: 0.00\nPhase: 0.00°', ...
             'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], ...
             'FontSize', 10, 'Tag', 'infoText', 'HorizontalAlignment', 'left');
end

function sliderCallback(src, ax, param)
    fig = ancestor(ax, 'figure');
    params = get(fig, 'UserData');
    
    % 获取滑块值并更新显示
    value = get(src, 'Value');
    valText = get(src, 'UserData');
    set(valText, 'String', sprintf('%.2f', value));
    
    % 更新参数
    switch param
        case 'I_RF'
            params.I_RF_phi = value;
        case 'Q_RF'
            params.Q_RF_phi = value;
        case 'PS5'
            params.PS5_fai = value;
    end
    
    % 更新绘图
    updatePlot(fig, ax, params);
end

function updatePlot(fig, ax, params)
    % 计算输出信号
    Eo = SEPP_IQ(params.IQ_IN, params.R, params.IQ_fai0, ...
                  params.up_fai0, params.I_RF_phi, params.PS2_fai, ...
                  params.down_fai0, params.Q_RF_phi, params.PS3_fai, ...
                  params.PS5_fai, params.k, params.b);
    
    % 提取实部和虚部
    I = real(Eo);
    Q = imag(Eo);
    
    % 保存参数
    set(fig, 'UserData', params);
    
    % 绘制星座图
    cla(ax);
    hold(ax, 'on');
    
    % 绘制单位圆
    theta = linspace(0, 2*pi, 100);
    plot(ax, cos(theta), sin(theta), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
    
    % 绘制星座点
    plot(ax, I, Q, 'o', 'MarkerSize', 12, 'LineWidth', 2, ...
        'MarkerEdgeColor', [0 0.447 0.741], 'MarkerFaceColor', [0 0.447 0.741]);
    plot(ax, [0 I], [0 Q], 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);
    
    hold(ax, 'off');
    
    % 更新信息显示
    infoText = findobj(fig, 'Tag', 'infoText');
    set(infoText, 'String', sprintf('I: %.2f\nQ: %.2f\nAmplitude: %.2f\nPhase: %.2f°', ...
        I, Q, abs(Eo), angle(Eo)*180/pi));
    
    % 保证图形为正方形
    axis(ax, 'equal');
    xlim(ax, [-1.2 1.2]);
    ylim(ax, [-1.2 1.2]);
    grid(ax, 'on');
end

%% 2to2 single-ended push-pull IQ
function Eo = SEPP_IQ(IQ_IN, R, IQ_fai0, ...
up_fai0, I_RF_phi, PS2_fai, ...
down_fai0, Q_RF_phi, PS3_fai, ...
PS5_fai, k, b)
    % 调用IQ函数
    [~,Eo] = IQ (IQ_IN,0,R, IQ_fai0, ...
    up_fai0,I_RF_phi,-I_RF_phi,PS2_fai,0, ... 
    down_fai0,Q_RF_phi,-Q_RF_phi,0,PS3_fai, ... 
    0,PS5_fai, ...
    k,b);
end

%% 2to2 IQ in layout
function Eo = layout_IQ(IQ_IN, R, IQ_fai0, ...
up_fai0, I_RF_phi, PS2_fai, ...
down_fai0, Q_RF_phi, PS3_fai, ...
PS5_fai, k, b)
    % 调用IQ函数
    [~,Eo] = IQ (IQ_IN,0,R, IQ_fai0, ...
    up_fai0,0,I_RF_phi,PS2_fai,0, ... 
    down_fai0,0,Q_RF_phi,0,PS3_fai, ... 
    0,PS5_fai, ...
    k,b);
end

%% 2to2 IQ 
function [Eo1,Eo2] = IQ (Ein1,Ein2,R, IQ_fai0, ...
up_fai0,up_CDM_fai1,up_CDM_fai2,up_PS_fai1,up_PS_fai2, ... 
down_fai0,down_CDM_fai1,down_CDM_fai2,down_PS_fai1,down_PS_fai2, ... 
IQ_up_PS_fai,IQ_down_PS_fai, ...
k,b)
    [E1,E2] = BS(R,Ein1,Ein2);
    [~,E3] = MZI(E1,0,R,up_fai0,up_CDM_fai1,up_CDM_fai2,up_PS_fai1,up_PS_fai2,k,b);
    [~,E4] = MZI(E2,0,R,down_fai0,down_CDM_fai1,down_CDM_fai2,down_PS_fai1,down_PS_fai2,k,b);
    E5 = E3*exp(1i*IQ_up_PS_fai)*exp(1i*IQ_fai0);
    E6 = E4*exp(1i*IQ_down_PS_fai);
    [Eo1,Eo2] = BS(R,E5,E6);
end

%% 2to2 MZI 
function [Eo1,Eo2] = MZI(Ein1,Ein2,R,fai0,CDM_fai1,CDM_fai2,PS_fai1,PS_fai2,k,b)
    % fai_0是上下臂固有相位差
    loss_factor1 = k*CDM_fai1+b;   % k is negative
    loss_factor2 = k*CDM_fai2+b;
    [E1,E2] = BS(R,Ein1,Ein2);
    E3 = E1*exp(1i*CDM_fai1)*sqrt(10^(-loss_factor1/10));
    E4 = E2*exp(1i*CDM_fai2)*sqrt(10^(-loss_factor2/10));
    E5 = E3*exp(1i*PS_fai1)*exp(1i*fai0);
    E6 = E4*exp(1i*PS_fai2);
    [Eo1,Eo2] = BS(R,E5,E6);
end

%% 2to2 BS (MMI)
function [Eo1,Eo2] = BS(R,Ein1,Ein2)
    T = 1-R;
    Eo1 = 1i*sqrt(R)*Ein1 + sqrt(T)*Ein2;
    Eo2 = sqrt(T)*Ein1 + 1i*sqrt(R)*Ein2;
end
