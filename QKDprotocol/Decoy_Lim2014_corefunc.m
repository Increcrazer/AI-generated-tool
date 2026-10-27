function [R_bitperpulse, R_bitpersecond, e_obs, phi_X, nX, detail] = Decoy_Lim2014_corefunc ... 
            (epsilon_cor, epsilon_sec, ...
             qX, N, f, gate_width, width_3dB, ...
             k, pk, ...
             L, alpha, eta_Bob_detect, ...
             e_mis_X, e_mis_Z, f_EC, ...
             p_ap, dc_count, deadtime)
    % 输入：
    % epsilon_cor, epsilon_sec: security criteria
    % qX: 选X基概率(efficient BB84), N: A发送的块长, f[Hz]: 频率, gate_width[s]: 门宽,
    % 3dBwidth[s]: 脉冲半高全宽
    % k = [0.7 0.1 0.0002]; 各个强度态的平均光子, pk = [0.5 0.25 0.25]; A端各个强度态的发送概率
    % L[km]: A到B距离, alpha[dB/km]: 光路损耗
    % eta_Bob_detect: Bob总传输/探测效率（四通道等效率），不含分基概率及软件窗。
    % e_mis_X(Z): A端态制备的X(Z)基错误率，用X(Z)基对比度反算得到, f_EC: 纠错效率
    % p_ap: 标量后脉冲概率，额外点击p_ap*D，其中一半计为误码；SNSPD当前设0。
    % dc_count[/s]: 四通道暗计数，顺序X0,X1,Z0,Z1，各基取单通道均值。
    % deadtime[s]: 标量有效死时间/屏蔽时间，每基共用连续计数率近似修正。
    % 探测计数模型与Weina/Decoy_Li2025_corefunc一致，AB选基概率均为qX。
    % 此简化模型不描述四通道独立恢复、时钟分箱或软件窗外原始负载。
    % 输出：
    % R_bitperpulse, R_bitpersecond:密钥率
    % e_obs: X基比特误码率（用于成码的是X基）
    % phi_X: 估计得到的用于成码的X基的相位误码率
    % nX: B接收到的X基块长
    
    validate_inputs(epsilon_cor,epsilon_sec,qX,N,f,gate_width,width_3dB, ...
        k,pk,L,alpha,eta_Bob_detect,e_mis_X,e_mis_Z,f_EC,p_ap,dc_count,deadtime);
    k = reshape(k,1,3); pk = reshape(pk,1,3);

    %% 定义光路到Bob解码的效率
    eta_ch = 10.^(-alpha.*L./10);
    [~, width_eff] = gaussian_pulse_analysis(width_3dB, gate_width); % 卡门宽带来的效率
    eta_sys = eta_ch.*eta_Bob_detect.*width_eff;
    [counts,detector] = calculate_counts(qX,N,f,gate_width,k,pk,eta_sys, ...
        width_eff,e_mis_X,e_mis_Z,p_ap,dc_count,deadtime);
    n_X = counts.n_X; n_Z = counts.n_Z;
    m_X = counts.m_X; m_Z = counts.m_Z;
    assert(all(m_X<=n_X) && all(m_Z<=n_Z) && sum(n_X+n_Z)<=N*(1+1e-12), ...
        'Lim2014:Counts','Need errors<=detections and total detections<=N.');
    nX = sum(n_X); nZ = sum(n_Z);
    mX = sum(m_X);
    mZ = sum(m_Z); 

    e_obs = mX./max(nX,realmin);

    %% 计算phi_X、SX_0、SX_1    
    nxk_plus = zeros(1,3);
    nxk_minus = zeros(1,3);
    nzk_plus = zeros(1,3);
    nzk_minus = zeros(1,3);
    mzk_plus = zeros(1,3);
    mzk_minus = zeros(1,3);

    for i = 1:3
       nxk_plus(i) =  calculate_nxk_plus(k(i), pk(i), n_X(i), nX, epsilon_sec);
       nxk_minus(i) = calculate_nxk_minus(k(i), pk(i), n_X(i), nX, epsilon_sec);
       nzk_plus(i) = calculate_nxk_plus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);
       nzk_minus(i) = calculate_nxk_minus(k(i), pk(i), n_Z(i), nZ, epsilon_sec);       
       mzk_plus(i) = calculate_mzk_plus(k(i), pk(i), m_Z(i), mZ, epsilon_sec);
       mzk_minus(i) = calculate_mzk_minus(k(i), pk(i), m_Z(i), mZ, epsilon_sec);
    end
    tau0=calculate_tau_n(k,pk,0); tau1=calculate_tau_n(k,pk,1);
    SX_0=calculate_SX_0(tau0,k(2),k(3),nxk_plus(2),nxk_minus(3));
    SX_0=max(0,min(nX,SX_0));
    SX_1=calculate_SX_1(tau0,tau1,k(1),k(2),k(3), ...
        nxk_plus(1),nxk_minus(2),nxk_plus(3),SX_0);
    SX_1=max(0,min(nX-SX_0,SX_1));

    SZ_0=calculate_SX_0(tau0,k(2),k(3),nzk_plus(2),nzk_minus(3));
    SZ_0=max(0,min(nZ,SZ_0));
    SZ_1=calculate_SX_1(tau0,tau1,k(1),k(2),k(3), ...
        nzk_plus(1),nzk_minus(2),nzk_plus(3),SZ_0);
    SZ_1=max(0,min(nZ-SZ_0,SZ_1));
             
    nuZ_1 = max(0,min(mZ,calculate_nu_Z_1(tau1, mzk_plus(2), mzk_minus(3), k(2), k(3))));
    phi_X = calculate_phi_X(SX_1, nuZ_1, SZ_1, epsilon_sec);
    
    %% 计算密钥率
    if isfinite(phi_X) && phi_X >= 0 && phi_X <= 0.5 && e_obs >= 0 && e_obs <= 0.5 % 添加有效性检查
        l = calculate_l(SX_0, SX_1, phi_X, f_EC.*binary_entropy(e_obs).*nX, epsilon_sec, epsilon_cor);
    else
        l = 0;
    end

    if ~isfinite(l) || ~isreal(l)
        l = 0;
    end
    l = floor(max(0,min(nX,l)));
    
    R_bitperpulse= l./N;
    R_bitpersecond = l./N.*f;
    detail = struct('l',l,'n_X',n_X,'n_Z',n_Z,'m_X',m_X,'m_Z',m_Z, ...
        'SX_0',SX_0,'SX_1',SX_1,'SZ_0',SZ_0,'SZ_1',SZ_1, ...
        'nuZ_1',nuZ_1,'phi_X',phi_X,'detector',detector);
end

function [counts,model]=calculate_counts(qX,N,f,w,k,pk,eta_sys,width_eff, ...
        eX,eZ,p_ap,dc_count,deadtime)
% Same simplified count equations as Weina/Decoy_Li2025_corefunc.
    dc_count=reshape(dc_count,1,4);
    p_dc=[mean(dc_count(1:2)),mean(dc_count(3:4))]*w;
    assert(all(p_dc<0.5),'Lim2014:DarkCounts', ...
        'Per-detector dark probability must be <0.5.');
    P=[qX^2,(1-qX)^2]; e=[eX,eZ];
    n=zeros(2,3); m=zeros(2,3); Cdt=zeros(1,2);
    primary=zeros(2,3); raw=zeros(2,3);
    signal_click=-expm1(-eta_sys*k);
    for basis=1:2
        D=P(basis)*pk.*(1-(1-2*p_dc(basis))*exp(-eta_sys*k));
        raw_click=D*(1+p_ap);
        Cdt(basis)=1/(1+f*sum(raw_click)*deadtime);
        Rk=Cdt(basis)*raw_click;
        errors=Cdt(basis)*(P(basis)*pk.*(p_dc(basis)+e(basis)*signal_click)+p_ap*D/2);
        assert(all(Rk<=P(basis)*pk),'Lim2014:Model', ...
            'Afterpulse model exceeds one detection per matched pulse.');
        n(basis,:)=N*Rk; m(basis,:)=N*errors;
        primary(basis,:)=D; raw(basis,:)=raw_click;
    end
    counts=struct('n_X',n(1,:),'n_Z',n(2,:),'m_X',m(1,:),'m_Z',m(2,:));
    model=struct('eta_sys',eta_sys,'width_eff',width_eff,'p_dc',p_dc, ...
        'Cdt',Cdt,'primary_click',primary,'raw_click',raw, ...
        'p_detect',n/N,'p_error',m/N,'model','Li2025_simplified_counts');
end

%%
function [sigma, energy_ratio] = gaussian_pulse_analysis(fwhm, gate_width)
% GAUSSIAN_PULSE_ANALYSIS 计算高斯脉冲的标准差和指定门宽内的能量占比
%
% 输入参数:
%   fwhm: 高斯脉冲的3dB全宽（FULL WIDTH AT HALF MAXIMUM）
%   gate_width: 以中心为基准的门的宽度（总宽度为 ±gate_width./2）
%
% 输出参数:
%   sigma: 高斯脉冲的标准差
%   energy_ratio: 在 gate_width 范围内的能量占比（0~1）

    % 1. 根据 FWHM 计算标准差 sigma
    % FWHM = 2 .* sqrt(2 .* log(2)) .* sigma ≈ 2.355 .* sigma
    sigma = fwhm ./ (2 .* sqrt(2 .* log(2)));  % 精确计算

    % 2. 计算门宽内的能量占比（使用误差函数 erf 计算积分）
    % 高斯函数在 [-a, a] 内的积分 = erf(a ./ (sqrt(2).*sigma))
    a = gate_width ./ 2;  % 半宽
    energy_ratio = erf(a ./ (sqrt(2) .* sigma));
end

%%
function l = calculate_l(SX_0, SX_1, phiX, lambda_EC, epsilon_sec, epsilon_cor)
    l = SX_0 + SX_1 - SX_1 .* binary_entropy(phiX) - lambda_EC - 6 .* log2(21./epsilon_sec) - log2(2./epsilon_cor);
end

function SX_0 = calculate_SX_0(tau_0, mu_2, mu_3, nX_mu2_plus, nX_mu3_minus)
    SX_0 = tau_0 .* (mu_2 .* nX_mu3_minus - mu_3 .* nX_mu2_plus) ./ (mu_2 - mu_3);
end

function SX_1 = calculate_SX_1(tau_0, tau_1, mu_1, mu_2, mu_3, n_plus_X_mu1, n_minus_X_mu2, n_plus_X_mu3, SX_0)
    SX_1 = tau_1 .* mu_1 .* ...
          (n_minus_X_mu2 - n_plus_X_mu3 - ...
           ((mu_2.^2 - mu_3.^2) ./ mu_1.^2) .* (n_plus_X_mu1 - SX_0 ./ tau_0)) ./ ...
          (mu_1 .* (mu_2 - mu_3) -  mu_2.^2 + mu_3.^2);
end

function phi_X = calculate_phi_X(SX_1, nuZ_1, SZ_1, epsilon_sec)
    if SX_1 <= 0 || SZ_1 <= 0 || nuZ_1 <= 0 || nuZ_1./SZ_1 >= 0.5
        phi_X = 0.5; % Conservative abort at singular/out-of-domain sampling bounds.
        return;
    end
    phi_X = (nuZ_1 ./ SZ_1) + calculate_gamma(epsilon_sec, nuZ_1./SZ_1, SZ_1, SX_1);
    phi_X = min(0.5,phi_X);
end

function nu_Z_1 = calculate_nu_Z_1(tau1, m_Z_mu2, m_Z_mu3, mu2, mu3)
    nu_Z_1 = tau1 .* (m_Z_mu2 - m_Z_mu3) ./ (mu2 - mu3);
end

%%
function nxk_plus = calculate_nxk_plus(k, pk, nx_k, nx, epsilon_sec)
    nxk_plus = (exp(k)./pk) .* (nx_k + sqrt(nx./2 .* log(21./epsilon_sec)));
end

function nxk_minus = calculate_nxk_minus(k, pk, nx_k, nx, epsilon_sec)
    nxk_minus = (exp(k)./pk) .* (nx_k - sqrt(nx./2 .* log(21./epsilon_sec)));
end

%%
function mzk_plus = calculate_mzk_plus(k, pk, mz_k, mz, epsilon_sec)
    mzk_plus = (exp(k)./pk) .* (mz_k + sqrt((mz./2) .* log(21./epsilon_sec)));
end

function mzk_minus = calculate_mzk_minus(k, pk, mz_k, mz, epsilon_sec)
    mzk_minus = (exp(k)./pk) .* (mz_k - sqrt((mz./2) .* log(21./epsilon_sec)));
end

%%
function gamma = calculate_gamma(a, b, c, d)
    log_arg = log(c+d)-log(c)-log(d)-log1p(-b)-log(b) ...
        +2*log(21)-2*log(a);
    eqa = (1./c+1./d).*(1-b).*b./log(2).^2.*log_arg;
    if eqa >= 0
        gamma = sqrt(eqa);
    else
        gamma = inf;
    end
end

function validate_inputs(ec,es,q,N,f,w,fwhm,k,pk,L,alpha,eta,ex,ez,fec,ap,dc,dt)
    validateattributes(ec,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(es,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(q,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(N,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(f,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(w,{'numeric'},{'real','finite','scalar','positive','<=',1/f});
    validateattributes(fwhm,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(k,{'numeric'},{'real','finite','vector','numel',3,'nonnegative'});
    validateattributes(pk,{'numeric'},{'real','finite','vector','numel',3,'positive'});
    assert(k(1)>k(2)+k(3) && k(2)>k(3),'Lim2014:IntensityOrder', ...
        'Require mu1>mu2+mu3 and mu2>mu3>=0.');
    assert(abs(sum(pk)-1)<1e-12,'Lim2014:Probabilities','pk must sum to one.');
    validateattributes(L,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(alpha,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(eta,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
    validateattributes(ex,{'numeric'},{'real','finite','scalar','>=',0,'<=',0.5});
    validateattributes(ez,{'numeric'},{'real','finite','scalar','>=',0,'<=',0.5});
    validateattributes(fec,{'numeric'},{'real','finite','scalar','>=',1});
    validateattributes(ap,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
    validateattributes(dc,{'numeric'},{'real','finite','vector','numel',4,'nonnegative'});
    validateattributes(dt,{'numeric'},{'real','finite','scalar','nonnegative'});
end

function tau_n = calculate_tau_n(K, Pk, n)
tau_n = 0;
    for i = 1:length(K)
        term = exp(-K(i)) .* K(i).^n .* Pk(i) ./ factorial(n);
        tau_n = tau_n + term;
    end
end

function entropy = binary_entropy(p)
    % Use the continuous limits h(0)=h(1)=0; invalid inputs remain NaN.
    entropy = NaN(size(p));
    entropy(p == 0 | p == 1) = 0;
    interior = p > 0 & p < 1;
    entropy(interior) = -p(interior) .* log2(p(interior)) ...
        - (1-p(interior)) .* log2(1-p(interior));
end
