classdef MetalLine_CapacitanceFunctions_CHANG
    methods (Static)
        % 单接地平面电容计算代码
        function C = calculate_MG_capacitance(t, h, W, epsilon)
            % 计算中间变量B
            B = 1 + t / h;

            % 计算中间变量p
            p = 2 * B^2 - 1 + sqrt((2 * B^2 - 1)^2 - 1);

            % 计算eta的表达式
            term1 =  pi * W / (2 * h);
            term2 = (p + 1) / (2 * sqrt(p)) * (1 + log(4 / (p - 1)));
            term3 = -2 * atanh(1 / sqrt(p));
            eta = sqrt(p) *(term1 + term2 + term3);

            % 计算Delta
            Delta = max(eta, p);

            % 计算Rb
            if W/h >= 5
                Rb = eta + (p + 1) / 2 * log(Delta);
            elseif W/h >= 1
                Rb = eta + (p + 1) / 2 * log(Delta);
                Rb = Rb - sqrt((Rb - 1)*(Rb - p)) + (p + 1) * atanh(sqrt((Rb-p)/(Rb-1))) ... 
                - 2*sqrt(p) * atanh(sqrt((Rb-p)/(p*(Rb-1)))) + (pi*W)/(2*h) * sqrt(p);
            else
                error('模型失效：W/h=%.3f < 1 (要求W/h≥1)', W/h);
            end
                   
            % 计算Ra的对数部分
            ln_Ra_term1 = -1 - pi * W / (2 * h);
            ln_Ra_term2 = -((p + 1) / sqrt(p)) * atanh(1 / sqrt(p));
            ln_Ra_term3 = -log((p - 1) / (4 * p));
            ln_Ra = ln_Ra_term1 + ln_Ra_term2 + ln_Ra_term3;
            Ra = exp(ln_Ra);

            % 计算电容C
            C = (2 * epsilon / pi) * log(2 * Rb / Ra);
        end
        % 双接地平面电容计算代码
        function C = calculate_GMG_capacitance(t, h, W, d, epsilon)
            % 计算中间变量α和γ
            alpha = (h + d + t) / h;
            gamma = d / h;

            % 计算中间变量q
            term_q = alpha^2 - gamma^2 - 1;
            q = 0.5 * (term_q + sqrt(term_q^2 - 4 * gamma^2));

            % 计算中间变量p
            p = q^2 / gamma^2;

            % 计算ln(R_A)
            term1_lnRA = -pi * W / (2 * h);
            term2_lnRA = -2 * alpha * atanh(sqrt((p + q) / (p * (1 + q))));
            term3_lnRA = 2 * gamma * atanh(1 / sqrt(p));
            term4_lnRA = log(4 * p / (p - 1));
            ln_RA = term1_lnRA + term2_lnRA + term3_lnRA + term4_lnRA;
            RA = exp(ln_RA);

            % 计算ln(R_B)
            term1_lnRB = pi * W / (2 * h) + 2 * alpha * atanh(sqrt((1 + q) / (p + q)));
            term2_lnRB = gamma * log((p - 1) / 4);
            term3_lnRB = -2 * atanh(1 / sqrt(p));
            ln_RB = (term1_lnRB + term2_lnRB + term3_lnRB)/gamma;
            RB = exp(ln_RB);

            % 计算单位长度电容C/s
            C_per_s = (2 / pi) * log(RB / RA);

            % 乘以介电常数得到最终电容值
            C = epsilon * C_per_s;
        end
    end
end

