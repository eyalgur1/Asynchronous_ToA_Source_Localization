function [funML, biasV, s, T, times] = WLS(array_to_alg)
%% Set parameters %%
N = array_to_alg.array_info.N;
n = array_to_alg.array_info.n; 
sensors = array_to_alg.sensors;
bias = array_to_alg.bias;
realization = array_to_alg.realization;
ML = array_to_alg.ML;
t = realization.t;  % ToA measurements

G = [2*sensors', -2*t', -ones(N,1)];
b = sum(sensors.^2)' - (t.^2)';
F_phase1 = [G'*G, -G'*b; -b'*G, b'*b];

%% Run SDP %%
tic
cvx_begin sdp quiet
variable X(n+3, n+3) symmetric
minimize trace(F_phase1*X)
subject to
X(n+2, n+3) == trace(X(1:n, 1:n)) - X(n+1, n+1);
X(n+3, n+3) == 1;
X == semidefinite(n + 3);
cvx_end

T = X(n+1, n+3);
B = 2*diag(t - T);
invR = inv(B*B);
invRb = invR*b;
GinvRb = G'*invRb;
F_phase2 = [G'*(invR*G), -GinvRb; -GinvRb', b'*invRb];

cvx_begin sdp quiet
variable X(n+3, n+3) symmetric
minimize trace(F_phase2*X)
subject to
X(n+2, n+3) == trace(X(1:n, 1:n)) - X(n+1, n+1);
X(n+3, n+3) == 1;
X == semidefinite(n + 3);
cvx_end
times = toc;

T = X(n+1, n+3);
s = X(n+3, 1:n)';
funML = ML(s, T, t);
biasV = bias(s, T);