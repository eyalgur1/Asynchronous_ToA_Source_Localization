function [funML, biasV, s, T, times] = SDP(array_to_alg)
%% Set parameters %%
N = array_to_alg.array_info.N;
n = array_to_alg.array_info.n; 
sensors = array_to_alg.sensors;
bias = array_to_alg.bias;
realization = array_to_alg.realization;
ML = array_to_alg.ML;
t = realization.t;  % ToA measurements
U = [eye(N) ones(N, 1)];
U2 = U'*U;
tU = t*U;
I = eye(n);
%%%%%%%%%%%%%%%%%%%%%%%%


%% Run SDP %%
tic
cvx_begin sdp quiet
variables s(n) z h(N + 1) 
variable H(N + 1, N + 1) symmetric
minimize trace(U2*H)-2*trace(tU*h)
subject to
for i=1:N
    H(i,i) == norm(sensors(:, i))^2 - 2*s'*sensors(:, i) + z 
end
[H h; h' 1] == semidefinite(N + 2);
[I s;s' z] == semidefinite(n + 1);
cvx_end
times = toc;

T = h(N + 1);
funML = ML(s, T, t);
biasV = bias(s, T);