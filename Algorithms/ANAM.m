function [funML, biasV, s, T, times] = ANAM(array_to_alg, init_to_alg)
%% Set parameters %%
sFDPG = init_to_alg.s_inner; 
rFDPG = init_to_alg.rFDPG; 
s = init_to_alg.s0; 
max_iter = init_to_alg.max_iter;
N = array_to_alg.array_info.N; 
n = array_to_alg.array_info.n; 
const_mat = array_to_alg.array_info.const_mat; 
e1 = array_to_alg.array_info.e1; 
ind_u = array_to_alg.array_info.ind_u; 
ind_v = array_to_alg.array_info.ind_v;
ones_n = ones(n, 1);
sensors = array_to_alg.sensors;
p = array_to_alg.p;
sensors_ind_v = sensors(:, ind_v);
sensors_ind_u = sensors(:, ind_u);
bias = array_to_alg.bias;
realization = array_to_alg.realization;
ML = array_to_alg.ML;
t = realization.t; sum_t = sum(t);  % ToA measurements

% initial points
zeta = [s zeros(n, N - 1)];  % dual initial point
T = (1/N)*sum(t - sqrt(sum((s - sensors).^2)));  % set intial T (for funML and biasV calculation)

% set arrays for fun values (funML) and bias (biasV) along the iterations
funML = zeros(max_iter + 1,1); biasV = funML; times = zeros(max_iter, 1);
funML(1) = ML(s, T, t);
biasV(1) = bias(s, T);

tot_iter = 0;  % total iteration counter
out_iter = 0;  % outer iteration counter
%%%%%%%%%%%%%%%%%%%%%%%%


%% Run ANAM %%
while tot_iter < max_iter
    out_iter = out_iter + 1;
    time_a_start = tic;

    %%%% Update u and v %%%%
    u = (s - sensors); v = u;
    norm_u = sqrt(sum(u.^2)); nz_norm = (norm_u > 0);
    if sum(nz_norm) < N  % if s = pi
        not_nz_norm = not(nz_norm); E1 = repmat(e1, [1, sum(not(nz_norm))]);
        u(:, nz_norm) = -u(:, nz_norm)./norm_u(nz_norm);
        v(:, nz_norm) = -u(:, nz_norm);
        u(:, not_nz_norm) = E1;
        v(:, not_nz_norm) = E1;
    else  % s is not any sensor
        u = -u./norm_u;
        v = -u;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%


    %%%% Update A, L and z %%%%
    u_ind_u = u(:, ind_u); v_ind_v = v(:, ind_v);
    sum_uv = sum(reshape(kron(u_ind_u, ones_n').*kron(reshape(v_ind_v, 1, []), ones_n), n, n, []), 3);
    A = const_mat + (1/N)*sum_uv; A = 0.5*(A + A');
    invL = (2*min(eig(A)))/N;

    z1 = sum(sum(sensors_ind_v.*v_ind_v).*u_ind_u, 2);
    z2 = sum(sum(sensors_ind_u.*u_ind_u).*v_ind_v, 2);
    z3 = sum(t.*u, 2);
    z = -p +(1/N)*(z1 + z2) - z3;

    %%%% FDPG %%%%
    % FDPG Initialization
    xi = zeta; alphj = 1;
    FDPG_iter = sFDPG + 2^floor(out_iter/rFDPG) - 1;
    times(tot_iter + 1) = toc(time_a_start);

    % FDPG Iterative steps
    for j = 1:FDPG_iter

        % Set itertion counter and times
        if tot_iter >= max_iter
            break
        end
        tot_iter = tot_iter + 1;
        time_b_start = tic;

        % Set previous iterations
        alphj_prev = alphj;
        zeta_prev = zeta;

        % Update dual variable zeta (n×N matrix form)
        grad = 0.5*(A\(z + sum(xi, 2))) - sensors;
        c = xi - invL*grad;
        zeta = sum_t*c./max(N*sqrt(sum(c.^2)), sum_t);

        % Update step-size alphaj and auxiliary variable xi (n×N matrix form)
        alphj = (1 + sqrt(1 + 4*alphj^2))/2;
        xi = zeta + ((alphj_prev-1)/alphj)*(zeta - zeta_prev);

        % Update primal variable s and time offset T
        s = 0.5*(A\(z + sum(zeta, 2)));
        T = (1/N)*sum(t - sqrt(sum((s - sensors).^2)));

        % Update times, ML function value and squared bias
        times(tot_iter) = toc(time_b_start) + times(tot_iter);
        funML(tot_iter + 1) = ML(s, T, t);
        biasV(tot_iter + 1) = bias(s, T);
    end
    %%%%%%%%%%%%%%
end
end