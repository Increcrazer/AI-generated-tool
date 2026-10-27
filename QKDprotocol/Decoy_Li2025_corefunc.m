function [R_bitperpulse, R_bitpersecond, e_obs, phi_X, nX, phi_Z, detail] = Decoy_Li2025_corefunc ...
            (epsilon_cor, epsilon_sec, ...
             qX, N, f, gate_width, width_3dB, ...
             k, pk, ...
             L, alpha, eta_Bob_detect, ...
             e_mis_X, e_mis_Z, f_EC, ...
             p_ap, dc_count, deadtime, options)
%DECOY_LI2025_COREFUNC Li et al., Nature 640 (2025), Eqs. (10)-(21).
% DOI: 10.1038/s41586-025-08739-z.
% 输入顺序沿用 Decoy_Lim2014_corefunc；每次处理一个距离/数据块。
% k=[mu,nu,w], pk=[p_signal,p_decoy,p_vacuum]，均为1x3。
% L[km], alpha[dB/km], f[Hz], gate_width/width_3dB/deadtime[s]。
% qX: Alice 选X基概率；Bob默认同值；dc_count[/s]:四通道，前二X后二Z。
% e_mis_X/Z:基矢失配误码率；f_EC:纠错效率；eta_Bob_detect:Bob总效率。
%
% 输出前五项沿用Lim命名。区别：Li仅信号态成钥，e_obs为X基信号态QBER；
% nX仍为所有强度X基筛后总计数。phi_X/Z为各自基的相位误码率上界。
% detail包含计数、单光子期望/实际计数界、两基QBER和安全参数。
%
% options（可选，第19个输入）:
%   intensity_delta: 光强相对区间半宽，标量或[signal,decoy,vacuum]，默认0。
%   k_bounds: 3x2光强硬界[lower,upper]，覆盖intensity_delta。
%   qX_Bob: Bob选X基概率，默认qX。
%   epsilon_hat,epsilon_PA,xi: 必须一起指定，满足式(21)<=epsilon_sec。
%   counts: 实测计数struct(n_X,n_Z,m_X,m_Z)，各1x3，强度顺序同k。
%           此时绕过光路模型；计数为整数；N/f表示当前块发送时长。
%   lambda_EC: 实际纠错泄露[bit]；默认两基信号态f_EC*n*h(QBER)之和。
%
% 无counts时用探测模型均值作为典型观测计数，允许非整数，用于规划仿真。
% 安全估计采用原文式(16): m1_real_upper=M_signal，不复制keyComm的扣背景
% 或诱骗态错误数转移分支。光强区间失效概率不包含在式(21)中。
% 模型p_ap增加p_ap*D个点击，其中一半为误码；每基共用有效死时间修正。

    %% 输入参数检查
    if nargin<19 || isempty(options), options=struct(); end
    validate_inputs(epsilon_cor,epsilon_sec,qX,N,f,k,pk,f_EC,options);
    k=reshape(k,1,3); pk=reshape(pk,1,3);

    %% 基矢选择概率、安全参数和光源界
    qX_Bob=get_option(options,'qX_Bob',qX);
    validateattributes(qX_Bob,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    qZ=1-qX; qZ_Bob=1-qX_Bob;
    P_basis=[qX*qX_Bob,qZ*qZ_Bob];
    security=calculate_security(epsilon_cor,epsilon_sec,options);
    [k_bounds,coeff]=calculate_source_bounds(k,options);

    %% 光路效率、探测计数和误码计数
    if isfield(options,'counts')
        counts=validate_observed_counts(options.counts);
        input_mode='observed_counts'; model=struct();
    else
        validate_model_inputs(gate_width,width_3dB,L,alpha,eta_Bob_detect, ...
            e_mis_X,e_mis_Z,p_ap,dc_count,deadtime);
        eta_ch=10^(-alpha*L/10);
        [~,width_eff]=gaussian_pulse_analysis(width_3dB,gate_width);
        eta_sys=eta_ch*eta_Bob_detect*width_eff;
        [counts,model]=calculate_counts(qX,N,f,gate_width,k,pk,eta_sys, ...
            width_eff,e_mis_X,e_mis_Z,p_ap,dc_count,deadtime,qX_Bob);
        input_mode='model_typical_counts';
    end
    n_X=counts.n_X; n_Z=counts.n_Z;
    m_X=counts.m_X; m_Z=counts.m_Z;
    n=[n_X;n_Z]; m=[m_X;m_Z];
    assert(all(m(:)<=n(:)) && sum(n(:))<=N*(1+1e-12), ...
        'Li2025:Counts','Need 0<=errors<=detections and sum(detections)<=N.');
    nX=sum(n_X); nZ=sum(n_Z);

    %% 信号态QBER和纠错泄露量
    e_signal=[calculate_ratio(m(1,1),n(1,1)),calculate_ratio(m(2,1),n(2,1))];
    e_obs=e_signal(1);
    if isfield(options,'lambda_EC')
        lambda_EC=options.lambda_EC;
        validateattributes(lambda_EC,{'numeric'},{'real','finite','scalar','nonnegative'});
    else
        lambda_EC=f_EC*sum(n(:,1)'.*binary_entropy(e_signal));
    end

    %% 观测计数到期望计数界：式(3)-(6)，xi为单侧失败概率
    expected_lower=zeros(2,3); expected_upper=zeros(2,3);
    for basis=1:2
        for intensity=1:3
            [expected_lower(basis,intensity),expected_upper(basis,intensity)]= ...
                calculate_expected_bounds(n(basis,intensity),security.xi,N);
        end
    end
    %% 信号态单光子期望数与实际数下界：式(12)-(15)
    expected_n1_lower=zeros(1,2); n1_real_lower=zeros(1,2);
    for basis=1:2
        probability=pk*P_basis(basis);
        expected_n1_lower(basis)=calculate_n1_expected( ...
            expected_lower(basis,:),expected_upper(basis,:),probability,coeff);
        expected_n1_lower(basis)=min(expected_n1_lower(basis),expected_upper(basis,1));
        % (13)-(15): expected -> actual; different from Lim2014's direct bound.
        n1_real_lower(basis)=min(n(basis,1), ...
            calculate_real_lower(expected_n1_lower(basis),security.xi));
    end
    %% 单光子误码数与相位误码率：式(16)-(19)
    % Actual single-photon errors <= observed total SIGNAL errors.
    m1_real_upper=m(:,1)';
    e1_bit_upper=ones(1,2);
    for basis=1:2
        if n1_real_lower(basis)>0
            e1_bit_upper(basis)=min(1,m1_real_upper(basis)/n1_real_lower(basis));
        end
    end
    % (17)-(19): Z bit errors -> X phase errors; X bit errors -> Z phase errors.
    phi_X=calculate_phi(e1_bit_upper(2),n1_real_lower(2),n1_real_lower(1),security.xi);
    phi_Z=calculate_phi(e1_bit_upper(1),n1_real_lower(1),n1_real_lower(2),security.xi);
    %% 密钥长度和成码率：式(20)
    % sqrt(2)*epsilon_PA*epsilon_hat in denominator (NOT 2*sqrt(2)).
    privacy_bits=n1_real_lower.*(1-binary_entropy([phi_X,phi_Z]));
    correction_bits=2*(1-log2(epsilon_cor)) + ...
        4*(-0.5-log2(security.epsilon_PA)-log2(security.epsilon_hat));
    l_raw=sum(privacy_bits)-lambda_EC-correction_bits;
    l=max(0,floor(min(sum(n(:,1)),l_raw)));
    R_bitperpulse=l/N; R_bitpersecond=R_bitperpulse*f;

    %% 中间统计量与状态输出
    detail=struct('l',l,'l_raw',l_raw,'input_mode',input_mode, ...
        'n_X',n(1,:),'n_Z',n(2,:),'m_X',m(1,:),'m_Z',m(2,:), ...
        'nX',nX,'nZ',nZ,'e_signal',e_signal,'P_basis',P_basis, ...
        'k_bounds',k_bounds,'coeff',coeff, ...
        'expected_counts_lower',expected_lower,'expected_counts_upper',expected_upper, ...
        'n1_expected_lower',expected_n1_lower,'n1_real_lower',n1_real_lower, ...
        'm1_real_upper',m1_real_upper,'e1_bit_upper',e1_bit_upper, ...
        'phi_X',phi_X,'phi_Z',phi_Z,'privacy_bits',privacy_bits, ...
        'lambda_EC',lambda_EC,'correction_bits',correction_bits, ...
        'security',security,'model',model,'source_condition',coeff.valid);
    if ~coeff.valid
        detail.status='source_condition_failed';
    elseif any(n1_real_lower<=0)
        detail.status='insufficient_single_photon_counts';
    elseif l==0
        detail.status='no_positive_key';
    else
        detail.status='ok';
    end
end

function security=calculate_security(epsilon_cor,epsilon_sec,options)
% (21): epsilon_sec=2*epsilon_hat+epsilon_PA+4*sqrt(5*xi).
% Default split: 1/4 for 2*hat, 1/4 for PA, 1/2 for statistical terms.
    names={'epsilon_hat','epsilon_PA','xi'}; supplied=isfield(options,names);
    assert(all(supplied) || ~any(supplied),'Li2025:Security', ...
        'Provide epsilon_hat,epsilon_PA,xi together.');
    if all(supplied)
        epsilon_hat=options.epsilon_hat; epsilon_PA=options.epsilon_PA; xi=options.xi;
    else
        epsilon_hat=epsilon_sec/8; epsilon_PA=epsilon_sec/4;
        xi=(epsilon_sec/8)^2/5;
    end
    validateattributes(epsilon_hat,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(epsilon_PA,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(xi,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    sec_used=2*epsilon_hat+epsilon_PA+4*sqrt(5*xi);
    assert(sec_used<=epsilon_sec*(1+1e-12),'Li2025:Security','Eq.(21) exceeds epsilon_sec.');
    security=struct('epsilon_cor',epsilon_cor,'epsilon_sec_target',epsilon_sec, ...
        'epsilon_sec_used',sec_used,'epsilon_hat',epsilon_hat,'epsilon_PA',epsilon_PA, ...
        'xi',xi,'epsilon_total',2*(sec_used+epsilon_cor));
end

function [bounds,c]=calculate_source_bounds(k,options)
% (10)-(11): conservative product bounds e^(-upper)*lower^n/n!,
% not simply the Poisson probability evaluated at lower intensity.
    if isfield(options,'k_bounds')
        bounds=options.k_bounds;
        validateattributes(bounds,{'numeric'},{'real','finite','size',[3,2],'nonnegative'});
    else
        delta=get_option(options,'intensity_delta',0);
        validateattributes(delta,{'numeric'},{'real','finite','vector','nonnegative','<',1});
        assert(isscalar(delta) || numel(delta)==3,'Li2025:IntensityBounds', ...
            'intensity_delta must be scalar or have three entries.');
        if isscalar(delta), delta=repmat(delta,1,3); end
        delta=reshape(delta,1,3);
        bounds=[k(:).*(1-delta(:)),k(:).*(1+delta(:))];
    end
    assert(all(bounds(:,1)<=k(:) & k(:)<=bounds(:,2)), ...
        'Li2025:IntensityBounds','k_bounds must enclose nominal k.');
    assert(bounds(1,1)>bounds(2,2) && bounds(2,1)>bounds(3,2), ...
        'Li2025:IntensityBounds','Intervals must not overlap: signal>decoy>vacuum.');
    % Rows vacuum/decoy/signal; columns photon numbers 0/1/2/3.
    lower=zeros(3,4); upper=zeros(3,4);
    for source=1:3
        interval=bounds(4-source,:);
        lower(source,1)=exp(-interval(2)); upper(source,1)=exp(-interval(1));
        for photons=1:3
            lower(source,photons+1)=exp(-interval(2))*interval(1)^photons/factorial(photons);
            upper(source,photons+1)=exp(-interval(1))*interval(2)^photons/factorial(photons);
        end
    end
    a=lower(1,:); b=upper(2,:); cc=lower(3,:);
    A=a(1)*cc(3)-cc(1)*a(3);
    B=a(1)*b(3)-b(1)*a(3);
    C=cc(3)*b(1)-b(3)*cc(1);
    denominator=A*(a(1)*b(2)-b(1)*a(2))-B*(a(1)*cc(2)-cc(1)*a(2));
    % Sufficient all-n condition for Eq.(12): g_n=B*c_nL-A*b_nU+C*a_nL>=0.
    % g_2=0. Dividing by (nuU)^n/n! yields a positive sum of exponentials
    % minus a constant, hence discrete convexity: g_3>=0 implies all n>=3.
    g3=B*cc(4)-A*b(4)+C*a(4);
    c=struct('lower',lower,'upper',upper,'A',A,'B',B,'C',C, ...
        'denominator',denominator,'g3',g3, ...
        'valid',A>0 && B>0 && C>0 && denominator>0 && g3>=0);
end

function bound=calculate_n1_expected(nL,nU,p,c)
% Eq.(12), p includes the matched-basis factor; order signal/decoy/vacuum.
    if ~c.valid, bound=0; return; end
    a0L=c.lower(1,1); c1L=c.lower(3,2);
    numerator=c.A*a0L*nL(2)/p(2)-c.B*a0L*nU(1)/p(1)-c.C*a0L*nU(3)/p(3);
    bound=max(0,p(1)*c1L*numerator/c.denominator);
end

function [lower,upper]=calculate_expected_bounds(x,xi,maximum)
% (3)-(6) in log space: x*log(x/mu)+mu-x=log(1/xi).
    t=-log(xi);
    if x==0, lower=0; upper=min(maximum,t); return; end
    lo=0; hi=x;
    for i=1:100
        mid=(lo+hi)/2;
        if poisson_deviance(x,mid)>t, lo=mid; else, hi=mid; end
    end
    lower=lo;
    lo=x; hi=x+sqrt(2*x*t)+2*t;
    while poisson_deviance(x,hi)<t, hi=2*hi; end
    for i=1:100
        mid=(lo+hi)/2;
        if poisson_deviance(x,mid)>t, hi=mid; else, lo=mid; end
    end
    upper=min(maximum,hi);
end

function d=poisson_deviance(x,mu)
    if mu==0, d=Inf; return; end
    r=(mu-x)/x;
    if abs(r)<1e-3
        % Stable expansion of r-log(1+r), needed for N~1e20.
        d=x*r*r*(1/2+r*(-1/3+r*(1/4+r*(-1/5+r*(1/6+r*(-1/7+r/8))))));
    else
        d=x*(r-log1p(r));
    end
end

function lower=calculate_real_lower(mu,xi)
% (13)-(15): mu*[delta+(1-delta)*log(1-delta)]=log(1/xi).
    t=-log(xi);
    if mu<=t, lower=0; return; end
    lo=0; hi=1;
    for i=1:100
        delta=(lo+hi)/2;
        if delta<1e-3
            h=delta^2*(1/2+delta*(1/6+delta*(1/12+delta*(1/20+delta/30))));
        else
            h=delta+(1-delta)*log1p(-delta);
        end
        if mu*h>t, hi=delta; else, lo=delta; end
    end
    lower=max(0,(1-hi)*mu);
end

function phi=calculate_phi(bit_error,c,d,xi)
% (17)-(19). phi=0.5 means zero certifiable privacy in this basis.
% Eq.(19) is singular at b=0; do not set gamma=0 there.
    phi=0.5;
    if c<=0 || d<=0 || bit_error<=0 || bit_error>=0.5, return; end
    b=bit_error;
    log_argument=log(c+d)-log(2*pi)-log(c)-log(d)-log1p(-b)-log(b)-2*log(xi);
    if log_argument<=0, return; end
    gamma=sqrt((1/c+1/d)*(1-b)*b*log_argument);
    phi=min(0.5,b+gamma);
end

function [counts,model]=calculate_counts(qA,N,f,gate_width,k,pk,eta_sys, ...
    width_eff,eX,eZ,p_ap,dc_count,deadtime,qB)
% Simplified detector model, shared in form with the Lim2014 implementation.
    dc_count=reshape(dc_count,1,4);
    p_dc=[mean(dc_count(1:2)),mean(dc_count(3:4))]*gate_width;
    assert(all(p_dc<0.5),'Li2025:DarkCounts','Per-detector dark probability must be <0.5.');
    P=[qA*qB,(1-qA)*(1-qB)]; e=[eX,eZ];
    n=zeros(2,3); m=zeros(2,3); Cdt=zeros(1,2);
    signal_click=-expm1(-eta_sys*k);
    for basis=1:2
        D=P(basis)*pk.*(1-(1-2*p_dc(basis))*exp(-eta_sys*k));
        raw_click=D*(1+p_ap);
        Cdt(basis)=1/(1+f*sum(raw_click)*deadtime);
        Rk=Cdt(basis)*raw_click;
        errors=Cdt(basis)*(P(basis)*pk.*(p_dc(basis)+e(basis)*signal_click)+p_ap*D/2);
        assert(all(Rk<=P(basis)*pk),'Li2025:Model', ...
            'Afterpulse model exceeds one detection per matched pulse.');
        n(basis,:)=N*Rk; m(basis,:)=N*errors;
    end
    counts=struct('n_X',n(1,:),'n_Z',n(2,:),'m_X',m(1,:),'m_Z',m(2,:));
    model=struct('eta_sys',eta_sys,'width_eff',width_eff,'p_dc',p_dc,'Cdt',Cdt);
end

function [sigma,energy_ratio]=gaussian_pulse_analysis(fwhm,gate_width)
% Gaussian pulse fraction within the centered software window.
    sigma=fwhm/(2*sqrt(2*log(2)));
    energy_ratio=erf(gate_width/(2*sqrt(2)*sigma));
end

function validate_inputs(epsilon_cor,epsilon_sec,qX,N,f,k,pk,f_EC,options)
    validateattributes(epsilon_cor,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(epsilon_sec,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(qX,{'numeric'},{'real','finite','scalar','>',0,'<',1});
    validateattributes(N,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(f,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(k,{'numeric'},{'real','finite','vector','numel',3,'nonnegative'});
    validateattributes(pk,{'numeric'},{'real','finite','vector','numel',3,'positive'});
    validateattributes(f_EC,{'numeric'},{'real','finite','scalar','>=',1});
    assert(abs(sum(pk)-1)<1e-12,'Li2025:Probabilities','sum(pk) must equal 1.');
    assert(k(1)>k(2) && k(2)>k(3),'Li2025:Intensities','Require mu>nu>w>=0.');
    assert(isstruct(options) && isscalar(options),'Li2025:Options','options must be a scalar struct.');
    allowed={'intensity_delta','k_bounds','qX_Bob','epsilon_hat','epsilon_PA', ...
        'xi','counts','lambda_EC'};
    assert(all(ismember(fieldnames(options),allowed)),'Li2025:Options','Unknown options field.');
end

function validate_model_inputs(gate_width,width_3dB,L,alpha,eta_Bob_detect, ...
    eX,eZ,p_ap,dc_count,deadtime)
    validateattributes(gate_width,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(width_3dB,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(L,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(alpha,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(eta_Bob_detect,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
    validateattributes(eX,{'numeric'},{'real','finite','scalar','>=',0,'<=',0.5});
    validateattributes(eZ,{'numeric'},{'real','finite','scalar','>=',0,'<=',0.5});
    validateattributes(p_ap,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
    validateattributes(deadtime,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(dc_count,{'numeric'},{'real','finite','vector','numel',4,'nonnegative'});
end

function counts=validate_observed_counts(counts)
    names={'n_X','n_Z','m_X','m_Z'};
    assert(isstruct(counts) && isscalar(counts) && all(isfield(counts,names)), ...
        'Li2025:Counts','counts needs n_X,n_Z,m_X,m_Z.');
    for i=1:numel(names)
        value=counts.(names{i});
        validateattributes(value,{'numeric'}, ...
            {'real','finite','vector','numel',3,'nonnegative','integer'});
        counts.(names{i})=reshape(value,1,3);
    end
end

function result=binary_entropy(p)
    result=zeros(size(p)); interior=p>0 & p<1;
    result(interior)=-p(interior).*log2(p(interior))-(1-p(interior)).*log2(1-p(interior));
end

function ratio=calculate_ratio(errors,counts)
    if counts==0, ratio=0; else, ratio=errors/counts; end
end

function value=get_option(options,name,default)
    if isfield(options,name), value=options.(name); else, value=default; end
end
