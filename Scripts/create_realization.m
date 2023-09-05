function [realization, realization_info] = create_realization(s_real, T_real, sensors, N, sigma)

noises = sigma*randn(1, N);  % measurement noises
t = sqrt(sum((s_real - sensors).^2)) + T_real + noises;  % ToA measurements

while_count = 0;
while sum(t <= 0) >= 1  % always positive measurements
    while_count = while_count + 1;
    noises = sigma*randn(1, N);
    t = sqrt(sum((s_real - sensors).^2)) + T_real + noises;
    if while_count == 10000
        error('The given hyper-parameters sigma or offset_factor do not generate positive ToA measurements')
    end
end

realization.t = t;
realization.noises = noises;

ML = @(s, T, t)0.5*sum((sqrt(sum((s - sensors).^2)) + T - t).^2);  % ML function value (depends on each realization because of t)
realization_info.ML = ML;