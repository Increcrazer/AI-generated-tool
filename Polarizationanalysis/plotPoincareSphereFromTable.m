function plotPoincareSphereFromTable(filename, varargin)
    % 读取数据并支持选择工作表
    % 用法: 
    % plotPoincareSphereFromTable(filename) - 读取第一个工作表
    % plotPoincareSphereFromTable(filename, 'Sheet', 'SheetName') - 读取指定工作表
    % plotPoincareSphereFromTable(filename, 'Sheet', 2) - 读取第二个工作表
    
    % 解析输入参数
    p = inputParser;
    addRequired(p, 'filename', @ischar);
    addParameter(p, 'Sheet', 1, @(x) isnumeric(x) || ischar(x));
    parse(p, filename, varargin{:});
    
    % 检查文件类型
    [~, ~, ext] = fileparts(filename);
    if strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
        % 读取Excel文件
        if ischar(p.Results.Sheet)
            data = readtable(filename, 'Sheet', p.Results.Sheet);
        else
            % 获取所有工作表名
            [~, sheets] = xlsfinfo(filename);
            if p.Results.Sheet <= length(sheets)
                data = readtable(filename, 'Sheet', sheets{p.Results.Sheet});
            else
                error('工作表索引超出范围');
            end
        end
    else
        % 读取其他格式文件(CSV等)
        data = readtable(filename);
    end
   
    voltage = data.Voltage_mV;

    % 提取S1, S2, S3数据
    S1 = data.S1;
    S2 = data.S2;
    S3 = data.S3;
    
    % 创建图形窗口
    fig = figure('Name', '邦加球上的Stokes参数轨迹', 'NumberTitle', 'off', ...
           'Color', 'white', 'Position', [100, 100, 800, 800]);
    hold on;
    axis equal;
    view(3);
    
    % 绘制更美观的邦加球
    [x, y, z] = sphere(100);
    h = surf(x, y, z, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
             'FaceColor', [0.6 0.8 1], 'FaceLighting', 'gouraud');
    
    % 添加半透明网格线
    for r = 0.2:0.2:0.8
        [x, y, z] = sphere(50);
        surf(r*x, r*y, r*z, 'FaceAlpha', 0.05, 'EdgeColor', 'none', ...
             'FaceColor', [0.7 0.7 0.7]);
    end
    
    % 绘制坐标轴和标签
    line([-1.2 1.2], [0 0], [0 0], 'Color', [0.8 0 0], 'LineWidth', 1.5); % S1轴
    line([0 0], [-1.2 1.2], [0 0], 'Color', [0 0.6 0], 'LineWidth', 1.5); % S2轴
    line([0 0], [0 0], [-1.2 1.2], 'Color', [0 0 0.8], 'LineWidth', 1.5); % S3轴
    
    % 添加坐标轴标签
    text(1.25, 0, 0, 'S1', 'Color', [0.8 0 0], 'FontSize', 14, 'FontWeight', 'bold');
    text(0, 1.25, 0, 'S2', 'Color', [0 0.6 0], 'FontSize', 14, 'FontWeight', 'bold');
    text(0, 0, 1.25, 'S3', 'Color', [0 0 0.8], 'FontSize', 14, 'FontWeight', 'bold');
    
    % 绘制Stokes参数点（不连线）
    scatter3(S1, S2, S3, 60, voltage, 'filled', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % 添加颜色条表示电压变化
    cb = colorbar;
    cb.Label.String = '电压 (mV)';
    cb.Label.FontSize = 12;
    caxis([min(voltage) max(voltage)]);
    
    % 设置图形标题和标签
    if ischar(p.Results.Sheet)
        titleStr = sprintf('邦加球上的Stokes参数分布 - %s', p.Results.Sheet);
    else
        titleStr = '邦加球上的Stokes参数分布';
    end
    title(titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('S1', 'FontSize', 12);
    ylabel('S2', 'FontSize', 12);
    zlabel('S3', 'FontSize', 12);
    
    % 添加光照效果
    light('Position',[1 1 1],'Style','infinite');
    light('Position',[-1 -1 -1],'Style','infinite');
    material dull;
    
    % 调整视角和轴属性
    view(135, 30);
    axis off;
    grid off;
    
    % 添加邦加球赤道线（可选）
    theta = linspace(0, 2*pi, 100);
    plot3(cos(theta), sin(theta), zeros(size(theta)), 'k:', 'LineWidth', 0.5);
    
    hold off;
end
