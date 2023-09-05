function [funML, biasV, sk, Tk, times] = PAMP(array_to_alg, init_to_alg)

s0 = init_to_alg.s0;
T0 = init_to_alg.T0;
max_iter = init_to_alg.max_iter;
N = array_to_alg.array_info.N;
n = array_to_alg.array_info.n;
sensors = array_to_alg.sensors;
bias = array_to_alg.bias;
realization = array_to_alg.realization;
ML = array_to_alg.ML;
t = realization.t;  % ToA measurements
sPAMP = init_to_alg.s_inner;
alpha_min = init_to_alg.alpha_min;
alpha_max = init_to_alg.alpha_max;

funML = zeros(max_iter + 1,1); biasV = funML; times = zeros(max_iter, 1);

T_init_for_vals = (1/N)*sum(t - sqrt(sum((s0 - sensors).^2)));


funML(1) = ML(s0, T_init_for_vals, t);
biasV(1) = bias(s0, T_init_for_vals);

sk = s0;
Tk = T0;

g0 = zeros(n, 1); % initialize the gradient
for i = 1:N
    gi_part = sensors(:, i) - s0;
    norm_gi_part = norm(gi_part);
    if norm_gi_part > 0
        gi = gi_part./norm_gi_part;
    else
        gi = 0;
    end
    g0 = g0 + (t(i) - T0 - norm(sensors(:, i) - s0))*gi;
end
gl = g0;


out_iter = 0;
tot_iter = 0;

while tot_iter < max_iter
    time_start = tic;

    Tk_part = 0;
    for i=1:N
        Tk_part = Tk_part + (t(i) - norm(sensors(:, i) - sk));
    end
    Tk = (1/(N+1))*(Tk_part + Tk);


    out_iter = out_iter + 1;


    alphal = alpha_min;
    Ul = 1;
    sl = sk;

    for l=1:sPAMP
        if tot_iter >= max_iter
            break
        end

        sl_prev = sl;  % first inner iteration is the last outer

        % Claculate current gradinet, and stop is ||grad||<=eps^k
        gl_prev = gl;
        gl = zeros(n, 1); % initialize the gradient
        for i = 1:N
            gi_part = sensors(:, i) - sl;
            norm_gi_part = norm(gi_part);
            if norm_gi_part > 0
                gi = gi_part./norm_gi_part;
            else
                gi = 0;
            end
            gl = gl + (t(i) - Tk - norm(sensors(:, i) - sl))*gi;
        end

        is_back = 0;
        if norm(gl) >= eps^out_iter
            zl = sl - alphal*gl;

            phi_zl = 0;
            phi_sl = 0;
            for i = 1:N
                phi_zl = phi_zl + 0.5*(t(i) - Tk - norm(sensors(:, i) - zl))^2;
                phi_sl = phi_sl + 0.5*(t(i) - Tk - norm(sensors(:, i) - sl))^2;
            end
            delta_phi = phi_zl - phi_sl + alphal*norm(gl)^2;
            p = exp(-delta_phi/Ul);
            r = exp(-1) + (exp(-1/1)-exp(-1)).*rand(1);

            if p >= r
                sl = zl;
            else

                step_armijo = 1;
                sln = sl - step_armijo*gl;

                phi_sln = 0;
                for i = 1:N
                    phi_sln =  phi_sln + 0.5*(t(i) - Tk - norm(sensors(:, i) - sln))^2;
                end

                
                while phi_sln > phi_sl - 0.5*step_armijo*norm(gl)^2
                    is_back = 1;
                    if tot_iter >= max_iter
                        break
                    end
                    tot_iter = tot_iter + 1;
                    step_armijo = 0.5*step_armijo;
                    sln = sl - step_armijo*gl;

                    phi_sln = 0;
                    for i = 1:N
                        phi_sln = phi_sln + 0.5*(t(i) - Tk - norm(sensors(:, i) - sln))^2;
                    end

                    times(tot_iter) = toc(time_start);
                    funML(tot_iter + 1) = ML(sln, Tk, t);
                    biasV(tot_iter + 1) = bias(sln, Tk);
                end
                sl = sln;

            end

            alphal_BB = (norm(sl - sl_prev)^2)/((sl - sl_prev)'*(gl - gl_prev));

            alphal = max(alpha_min, min(alphal_BB, alpha_max));
            Ul = 0.5*Ul;
        end

        if is_back == 0
            tot_iter = tot_iter + 1;
            times(tot_iter) = toc(time_start);
            funML(tot_iter + 1) = ML(sl, Tk, t);
            biasV(tot_iter + 1) = bias(sl, Tk);
        end
    end

    sk = sl;
end
