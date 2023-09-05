clear; close all; cvx_solver SeDuMi; format compact

%% Hyper-parameters %%
%%%%%%%%%%%%%%%%%%%%%%

rng(8)
N = 15; % number of sensors
n = 2;  % dimensionarra
offset_factor = 0.04;  % time offset in Unifrom[-offset_factor,offset_factor]
sigma = logspace(-4,-1,6);  % noise factor for ToA measurements
s_inner_array = [15, 30];  % FDPG inner iterations (s)
rFDPG = 1000;  % FDPG inner iterations (r)
max_iter = 20000;  % total number of iterations (including inner) of ANAM
num_R = 50;  % number of measurement realizations for each array
num_A = 50;  % number of random arrays in [0,1]^n
start_point_convex_hull = 0;
% methods_to_run = {'ANAM','FPI','PAMP','SDP','WLS'};  % Required methods to run. The correct order is: ANAM -> FPI -> PAMP -> SDP (any can be omitted as long as the order is kept correct)
methods_to_run = {'ANAM','FPI','PAMP'}; 
alpha_min = 0.001;
alpha_max = 1;


%% Initialization %%
%%%%%%%%%%%%%%%%%%%%

len_s = length(s_inner_array); len_m = length(methods_to_run);
if sum(contains(methods_to_run, 'ANAM')) == 1
    if sum(contains(methods_to_run, 'PAMP')) == 1
        len_met = 2*len_s + len_m - 2;
    else
        len_met = len_s + len_m - 1;
    end
else
    len_met = len_m;
end
names_to_legend = cell(len_met, 1);
met_count = 0;

% set structure of num_A arrays with num_R realizations for each array
[arrays, arrays_info] = create_arrays(N, n, offset_factor, sigma, num_A, num_R);
out_struct = struct;


%% Run Algorithms %%
%%%%%%%%%%%%%%%%%%%%

for m = 1:len_m
    met = methods_to_run{m};

    switch met

        case 'ANAM'
            for j = 1:len_s
                met_count = met_count + 1; s_inner = s_inner_array(j); names_to_legend{met_count} = ['$\textrm{A-NAM},\ s=', num2str(s_inner),'$']; name_met = ['ANAM_s', num2str(s_inner)]; run_script_alg

                for arr = 1:num_A
                    run_script_array

                    for ss=1:length(sigma)
                        run_script_sig

                        for r = 1:num_R
                            fprintf(['ANAM s=',num2str(s_inner), ' || array=', num2str(arr), ' || sigma=', num2str(sig), ' || realization=', num2str(r), '\n']); run_script_real
                            [funML, biasV, s, T, times] = ANAM(array_to_alg, init_to_alg); run_script_save
                        end
                    end
                end
            end


        case 'FPI'
            met_count = met_count + 1; names_to_legend{met_count} = '$\textrm{FPI}$'; name_met = 'FPI'; run_script_alg

            for arr = 1:num_A
                run_script_array

                for ss=1:length(sigma)
                    run_script_sig

                    for r = 1:num_R
                        fprintf(['FP || array=', num2str(arr), ' || sigma=', num2str(sig), ' || realization=', num2str(r), '\n']); run_script_real
                        [funML, biasV, s, T, times] = FPI(array_to_alg, init_to_alg); run_script_save
                    end
                end
            end


        case 'PAMP'
            for j = 1:len_s
                met_count = met_count + 1; s_inner = s_inner_array(j); names_to_legend{met_count} = ['$\textrm{PAMP},\ s=', num2str(s_inner),'$']; name_met = ['PAMP_s', num2str(s_inner)]; run_script_alg
                out_struct.(name_met).alpha_min = alpha_min; out_struct.(name_met).alpha_max = alpha_max;

                for arr = 1:num_A
                    run_script_array

                    for ss=1:length(sigma)
                        run_script_sig

                        for r = 1:num_R
                            fprintf(['PAMP s=',num2str(s_inner), ' || array=', num2str(arr), ' || sigma=', num2str(sig), ' || realization=', num2str(r), '\n']); run_script_real
                            [funML, biasV, s, T, times] = PAMP(array_to_alg, init_to_alg); run_script_save
                        end
                    end
                end
            end


        case 'SDP'
            met_count = met_count + 1; names_to_legend{met_count} = '$\textrm{SDP}$'; name_met = 'SDP'; run_script_alg

            for arr = 1:num_A
                run_script_array

                for ss=1:length(sigma)
                    run_script_sig

                    for r = 1:num_R
                        fprintf(['SDP || array=', num2str(arr), ' || sigma=', num2str(sig), ' || realization=', num2str(r), '\n']); run_script_real
                        [funML, biasV, s, T, times] = SDP(array_to_alg); run_script_save
                    end
                end
            end


        case 'WLS'
            met_count = met_count + 1; names_to_legend{met_count} = '$\textrm{WLS}$'; name_met = 'WLS'; run_script_alg

            for arr = 1:num_A
                run_script_array

                for ss=1:length(sigma)
                    run_script_sig

                    for r = 1:num_R
                        fprintf(['WLS || array=', num2str(arr), ' || sigma=', num2str(sig), ' || realization=', num2str(r), '\n']); run_script_real
                        [funML, biasV, s, T, times] = WLS(array_to_alg); run_script_save
                    end
                end
            end
    end
end


%% Generate Plots

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1)
set(groot,'defaultAxesFontSize',14)

folderName = 'output/plots';
if ~exist(folderName, 'dir')
    mkdir(folderName)
end

formattedDateTime = string(char(datetime('now','Format','yyyy_MM_dd_hh_mm_ss')));

% ML function value and run time table (table for sig == 1)
out_struct = plot_ML_fun(sigma, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, max_iter, names_to_legend, N, formattedDateTime);

% Bias
plot_Bias(sigma, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, max_iter, names_to_legend, N, formattedDateTime)

% Plot RMSE and CRLB
plot_RMSE_CRLB(n, sigma, len_met, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, arrays_info, names_to_legend, N, formattedDateTime);

arrays.out_struct = out_struct;
save("output/output_TOA_"+string(N)+"_"+formattedDateTime+".mat", 'arrays', '-v7.3')