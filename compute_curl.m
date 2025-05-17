function curl_vec = compute_curl(Fx, Fy, Fz, x, y, z)
% 计算向量场 F = [Fx, Fy, Fz] 的旋度
% 输入:
%   Fx, Fy, Fz - 向量场的x、y、z分量（关于x,y,z的符号表达式）
%   x, y, z    - 符号变量
% 输出:
%   curl_vec - 旋度向量 [curl_x; curl_y; curl_z]

% 计算旋度的各分量
curl_x = diff(Fz, y) - diff(Fy, z);
curl_y = diff(Fx, z) - diff(Fz, x);
curl_z = diff(Fy, x) - diff(Fx, y);

% 组合成旋度向量
curl_vec = [curl_x; curl_y; curl_z];
end
