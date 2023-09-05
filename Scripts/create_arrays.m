function [arrays, arrays_info] = create_arrays(N, n, offset_factor, sigma, num_A, num_R)

arrays = struct; arrays_info = struct;
I = eye(n); e1 = I(:,1);

arrays_info.general.N = N;
arrays_info.general.n = n;
arrays_info.general.const_mat = ((N-1)/2)*I;
arrays_info.general.e1 = e1;
arrays.general.offset_factor = offset_factor;

% indices for pairs of u and v such that i<j
ind_uv = nchoosek(1:N, 2);
ind_u = ind_uv(:, 1)';
ind_v = ind_uv(:, 2)';
arrays_info.general.ind_u = ind_u;
arrays_info.general.ind_v = ind_v;

Jac = @(s,sensors)[((s - sensors)./sqrt(sum((s - sensors).^2)))' ones(N, 1)];  % Jacobian matrix for CRLB

for arr = 1:num_A   % generate num_A random arrays
    s_real = rand(n, 1);  % real source location
    arrays.(['array',num2str(arr)]).s_real = s_real;

    T_real = -offset_factor + 2*offset_factor*rand(1);  % real time-offset
    arrays.(['array', num2str(arr)]).T_real = T_real;  % real time-offset

    sensors = rand(n, N);  % sensors
    arrays.(['array', num2str(arr)]).sensors = sensors;
%     if n == 3
%         figure(1000+arr); hold on
%         scatter3(s_real(1), s_real(2), s_real(3), 'red')
%         scatter3(sensors(1, :), sensors(2, :),sensors(3, :),'blue')
%         hold off
%     elseif n==2
%         figure(1000+arr); hold on
%         scatter(s_real(1), s_real(2), 'red')
%         scatter(sensors(1, :), sensors(2, :), 'blue')
%         hold off
%     end

    arrays_info.(['array', num2str(arr)]).bias = @(s, T)norm(s - s_real)^2 + norm(T - T_real)^2;  % squared norm of bias function (each for every array)
    arrays_info.(['array', num2str(arr)]).p = (1/N - 1)*sum(sensors, 2);

    J = Jac(s_real,sensors);

    % generate the realizations for each array
    for ss = 1:length(sigma)
        sig = sigma(ss);
        reals = cell(num_R,1); reals_info = reals;  % empty cell for realizations (num_R realizations for each array)
        for r = 1:num_R
            [reals{r}, reals_info{r}] = create_realization(s_real, T_real, sensors, N, sig);
        end
        arrays.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations = reals;
        arrays_info.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations = reals_info;

        FIM = ((1/sig)^2)*(J'*J);
        arrays_info.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).CRLB =  sqrt(trace(inv(FIM)));

    end

end