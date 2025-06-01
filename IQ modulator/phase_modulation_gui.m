function phase_modulation_gui
    % 调制器由于制造工艺导致I、Q、IQ都具有不同初始相位
    % 在不同I、Q、IQ初始相位下扫描PS值，看输出光功率情况
    % 创建主窗口
    fig = figure('Name', '相位调制分析', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 1000, 700]);
    
    % 默认参数（存储在figure的UserData中）
    gui_data.IQ_IN = 1;
    gui_data.R = 0.5;
    gui_data.I_RF_phi = 0;
    gui_data.Q_RF_phi = 0;
    gui_data.k = -1.5;
    gui_data.b = 0;
    set(fig, 'UserData', gui_data);
    
    % 创建滑块和标签（调整了垂直位置）
    create_slider(fig, 'IQ_fai0', 'IQ相位 (IQ_fai0)', 0, 2*pi, 0, 0.01, [0.1, 0.80, 0.8, 0.03]);
    create_slider(fig, 'up_fai0', '上臂相位 (up_fai0)', 0, 2*pi, 0, 0.01, [0.1, 0.70, 0.8, 0.03]);
    create_slider(fig, 'down_fai0', '下臂相位 (down_fai0)', 0, 2*pi, 0, 0.01, [0.1, 0.60, 0.8, 0.03]);
    
    % 创建绘图区域（调整了位置）
    ax = axes('Parent', fig, 'Position', [0.1, 0.15, 0.85, 0.4]);
    
    % 初始绘图
    update_plot(fig);
    
    % 添加标题
    uicontrol('Style', 'text', 'String', '相位调制分析器', ...
              'Position', [350, 630, 300, 30], 'FontSize', 16, ...
              'FontWeight', 'bold', 'BackgroundColor', get(fig, 'Color'));
end

function create_slider(fig, name, label, min_val, max_val, init_val, step, position)
    % 创建标签
    uicontrol('Style', 'text', 'String', label, ...
              'Position', [position(1)*1000-100, position(2)*700+15, 200, 20], ...
              'FontSize', 10, 'HorizontalAlignment', 'right', ...
              'BackgroundColor', get(fig, 'Color'));
    
    % 创建滑块
    slider = uicontrol('Style', 'slider', 'Min', min_val, 'Max', max_val, ...
                       'Value', init_val, 'SliderStep', [step, step*10], ...
                       'Position', position.*[1000, 700, 1000, 700], ...
                       'Tag', name, 'Callback', {@slider_callback, fig});
    
    % 创建数值显示
    uicontrol('Style', 'text', 'String', sprintf('%.2f', init_val), ...
              'Position', [position(1)*1000+position(3)*1000+10, position(2)*700, 60, 20], ...
              'Tag', [name '_value'], 'BackgroundColor', get(fig, 'Color'));
end

function slider_callback(src, ~, fig)
    % 更新显示
    tag = get(src, 'Tag');
    value = get(src, 'Value');
    set(findobj(fig, 'Tag', [tag '_value']), 'String', sprintf('%.2f', value));
    
    % 更新绘图
    update_plot(fig);
end

function update_plot(fig)
    % 获取存储的参数
    gui_data = get(fig, 'UserData');
    
    % 获取当前滑块值
    IQ_fai0 = get(findobj(fig, 'Tag', 'IQ_fai0'), 'Value');
    up_fai0 = get(findobj(fig, 'Tag', 'up_fai0'), 'Value');
    down_fai0 = get(findobj(fig, 'Tag', 'down_fai0'), 'Value');
    
    % 创建相位范围
    phase_points = linspace(0, 2*pi, 100);
    
    % 预分配结果数组
    I_PS2 = zeros(size(phase_points));
    I_PS3 = zeros(size(phase_points));
    I_PS5 = zeros(size(phase_points));
    
    % 计算每个相位点的强度
    for i = 1:length(phase_points)
        % PS2_fai变化
        I_PS2(i) = abs(layout_IQ(gui_data.IQ_IN, gui_data.R, IQ_fai0, ...
                        up_fai0, gui_data.I_RF_phi, phase_points(i), ...
                        down_fai0, gui_data.Q_RF_phi, 0, ...
                        0, gui_data.k, gui_data.b))^2;
        
        % PS3_fai变化
        I_PS3(i) = abs(layout_IQ(gui_data.IQ_IN, gui_data.R, IQ_fai0, ...
                        up_fai0, gui_data.I_RF_phi, 0, ...
                        down_fai0, gui_data.Q_RF_phi, phase_points(i), ...
                        0, gui_data.k, gui_data.b))^2;
        
        % PS5_fai变化
        I_PS5(i) = abs(layout_IQ(gui_data.IQ_IN, gui_data.R, IQ_fai0, ...
                        up_fai0, gui_data.I_RF_phi, 0, ...
                        down_fai0, gui_data.Q_RF_phi, 0, ...
                        phase_points(i), gui_data.k, gui_data.b))^2;
    end
    
    % 获取坐标轴并绘制
    ax = findobj(fig, 'Type', 'axes');
    cla(ax);
    hold(ax, 'on');
    plot(ax, phase_points, I_PS2, 'r-', 'LineWidth', 2, 'DisplayName', 'I vs PS2\_fai');
    plot(ax, phase_points, I_PS3, 'g--', 'LineWidth', 2, 'DisplayName', 'I vs PS3\_fai');
    plot(ax, phase_points, I_PS5, 'b-.', 'LineWidth', 2, 'DisplayName', 'I vs PS5\_fai');
    hold(ax, 'off');
    
    % 设置图形属性
    xlabel(ax, '相位 (rad)');
    ylabel(ax, '强度 I');
    title(ax, '不同相位调制下的输出强度');
    legend(ax, 'show', 'Location', 'best');
    grid(ax, 'on');
    xlim(ax, [0, 2*pi]);
    ylim(ax, [0, max([I_PS2, I_PS3, I_PS5])*1.1]);
    
    % 添加当前相位值注释
    annotation_text = sprintf('当前相位值:\nIQ fai0 = %.2f\nup fai0 = %.2f\ndown fai0 = %.2f', ...
                             IQ_fai0, up_fai0, down_fai0);
    if isempty(findobj(fig, 'Type', 'annotation'))
        annotation(fig, 'textbox', [0.85, 0.85, 0.2, 0.15], 'String', annotation_text, ...
                   'FitBoxToText', 'on', 'BackgroundColor', 'white');
    else
        set(findobj(fig, 'Type', 'annotation'), 'String', annotation_text);
    end
end

%% 以下是原有的函数定义（保持不变）
%% 2to2 IQ in layout
function Eo = layout_IQ(IQ_IN, R, IQ_fai0, ...
    up_fai0, I_RF_phi, PS2_fai, ...
    down_fai0, Q_RF_phi, PS3_fai, ...
    PS5_fai, k, b)
    [~,Eo] = IQ(IQ_IN,0,R, IQ_fai0, ...
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

%% 2to2 BS
function [Eo1,Eo2] = BS(R,Ein1,Ein2)
    T = 1-R;
    Eo1 = 1i*sqrt(R)*Ein1 + sqrt(T)*Ein2;
    Eo2 =  sqrt(T)*Ein1 + 1i*sqrt(R)*Ein2;
end
