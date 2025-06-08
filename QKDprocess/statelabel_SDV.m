%% user set
filename = 'SDV-IM.csv';
bin_width = 16;
bin_number = 100000;
freq = 1.25 *10^9;  % [Hz]
pseudo_length = 32;
count_resol = 50;  % count resolution
arrange_list = char("SSDDSSVSVSSD");

%% find states of the pulses 
period = 1 /freq *10^12;    % [ps]
pseudo_time = 3*pseudo_length *period;   % [ps] 取3倍伪随机数长度

time = csvread(filename, 0, 0,[0,0,0,pseudo_time /bin_width]);
data = csvread(filename, 1, 0,[1,0,1,pseudo_time /bin_width]);
[~,index_first] = max(data(1:5*period /bin_width)); % 取5个脉冲内的第一个峰值
index_last = index_first + (pseudo_time -5*period) /bin_width;
index_list = index_first:period /bin_width:index_last;

pulse = zeros(1,length(index_list));
for i = 1:length(index_list)
    pulse(i) = sum(data(index_list(i)-period /bin_width /2:1:index_list(i)+period /bin_width /2));
end

log_pulse = log10(pulse);
figure(1);
[y,x] = hist(log_pulse,count_resol);   % y：每个柱子中包含的数据点数量（频数） x：每个柱子的中心位置坐标
bar(x,y);
xlabel('counts/log','FontName','Times New Roman','fontsize',18);
ylabel('probability','FontName','Times New Roman','fontsize',18);
title('Histogram');
nonzero_index = find(y); % find 函数用于查找数组中非零元素的索引

% find continuous numbers in nonzero_index and save in arrset
c1 = 1;
arrset = cell(0,0);
while (c1 < numel(nonzero_index))
    c2 = 0;
    while (c1+c2+1 <= numel(nonzero_index) && nonzero_index(c1)+c2+1 == nonzero_index(c1+c2+1))
        c2 = c2+1;
    end
    if (c1 >= 1)
        arrset= [arrset;(nonzero_index(c1:1:c1+c2))];
    end
    c1 = c1+c2+1;
end

if (numel(arrset) ~= 3)
    disp('请调整分辨率');
end

state1_range = [x(arrset{1}(1)),x(arrset{1}(end))];
state2_range = [x(arrset{2}(1)),x(arrset{2}(end))];
state3_range = [x(arrset{3}(1)),x(arrset{3}(end))];

% code the pulses
code_pulse = zeros(1,length(index_list));
for i = 1:length(index_list)
    if (log_pulse(i) >= state1_range(1) && log_pulse(i) <= state1_range(2))
        code_pulse(i) = 1;
    elseif (log_pulse(i) >= state2_range(1) && log_pulse(i) <= state2_range(2))
        code_pulse(i) = 2;
    elseif (log_pulse(i) >= state3_range(1) && log_pulse(i) <= state3_range(2))
        code_pulse(i) = 3;
    else
        code_pulse(i) = 0;
    end
end

% 实际上，SDV的情况毋庸置疑是dic6，不是说明有问题
% find pulse position
dic1 = struct('S',1,'D',2,'V',3);
dic2 = struct('S',1,'D',3,'V',2);
dic3 = struct('S',2,'D',1,'V',3);
dic4 = struct('S',2,'D',3,'V',1);
dic5 = struct('S',3,'D',1,'V',2);
dic6 = struct('S',3,'D',2,'V',1);

m = [];
for i = 1:6
    for j = 1:length(arrange_list)
        m = [m,eval(strcat('dic',num2str(i))).(arrange_list(j))];
    end
    eval([strcat('a',num2str(i)), '= m;']);
    m = [];
end

for i = 1:6
    m = eval(strcat('a',num2str(i)));
    for j = 1:length(code_pulse)-length(arrange_list)
        pulse_slice = code_pulse(j:j+length(arrange_list)-1);
        pulse_slice_index = find(code_pulse(j:j+length(arrange_list)-1));
        if (m(pulse_slice_index) == pulse_slice(pulse_slice_index))
            disp('Find position!');
            disp(strcat('Dictionary is dic',num2str(i)));
            disp(strcat('a(1) corresponds to code_pulse(',num2str(j),')'));
            arrange_list_order = j;
            break;
        else
            continue;
        end
    end
end


%%  plot and label state    
figure(2);
time = time(index_first-period /bin_width /2:index_last-period /bin_width /2);
data = data(index_first-period /bin_width /2:index_last-period /bin_width /2);
H = semilogy(time,data);
grid on;
xmin = time((arrange_list_order-1) *period /bin_width);
xmax = time((arrange_list_order+length(arrange_list)-1) *period /bin_width);
xlim([xmin, xmax]);
% ylim([0,5]);
xlabel('time/ps','FontName','Times New Roman','fontsize',18);
ylabel('counts/log','FontName','Times New Roman','fontsize',18);
title('SDV');
xpos = xmin+period/2:period:xmax-period/2;
ypos = 8000;
for i = 1:length(xpos)
    text(xpos(i),ypos,arrange_list(i));
end

