function [funML, biasV, s, T, times] = FPI(array_to_alg, init_to_alg)
%% Set parameters %%
s = init_to_alg.s0;
max_iter = init_to_alg.max_iter;
N = array_to_alg.array_info.N;
sensors = array_to_alg.sensors;
bias = array_to_alg.bias;
realization = array_to_alg.realization;
ML = array_to_alg.ML;
t = realization.t;  % ToA measurements

% initial points
diff = s - sensors;
T = (1/N)*sum(t - sqrt(sum(diff.^2)));

% set arrays for fun values (funML) and bias (biasV) along the iterations
funML = zeros(max_iter + 1,1); biasV = funML; times = zeros(max_iter, 1);
funML(1) = ML(s, T, t);
biasV(1) = bias(s, T);

iter = 0;  % total iteration counter
%%%%%%%%%%%%%%%%%%%%%%%%


%% Run FPI
while iter < max_iter
    time_start = tic;
    iter = iter + 1;

    % Calculation of FPI terms
    norm_diff = sqrt(sum(diff.^2));
    diff_normal = diff./norm_diff;
    part1 = (1 - 1/N)*sum((t.*diff_normal + sensors), 2);
    part2 = (1/N)*sum((t - norm_diff).*(sum(diff_normal, 2) - diff_normal), 2);

    % Update variable s and time offset T
    s = (1/(N-1))*(part1 + part2);
    diff = s - sensors;
    T = (1/N)*sum(t - sqrt(sum(diff.^2)));

    % Update times, ML function value and squared bias
    times(iter) = toc(time_start) + times(iter);
    funML(iter + 1) = ML(s, T, t);
    biasV(iter + 1) = bias(s, T);
end
