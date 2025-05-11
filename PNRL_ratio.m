function [PNRL_ratio] = PNRL_ratio(PN_ratio, RL_ratio)
    % PN_ratio RL_ratio为根据Bloch球夹角计算得出
    % 根据两态对比度算四态对比度
    PN_ER = 10^(-PN_ratio/10); % 为PN两态的平均误码率
    RL_ER = 10^(-RL_ratio/10); % 为RL两态的平均误码率
    PNRL_ER = (PN_ER + RL_ER)/2;
    PNRL_ratio = -10*log10(PNRL_ER);
end
