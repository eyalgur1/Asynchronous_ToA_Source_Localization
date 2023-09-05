function out_struct = plot_ML_fun(sigma, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, max_iter, names_to_legend, N, formattedDateTime)

total_run_time = [];
run_time_per_iter = [];

for ss = 1:length(sigma)
    sig = sigma(ss);
    figure(ss); hold on
    for m = 1:len_m
        met = methods_to_run{m};
        switch met

            case 'ANAM'
                for j = 1:len_s
                    avg_ML = zeros(max_iter + 1, 1);
                    if ss == 1
                        avg_times = zeros(max_iter, 1);
                    end
                    for arr = 1:num_A
                        for r = 1:num_R
                            avg_ML = avg_ML + (1/num_R)*out_struct.(['ANAM_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.fun_ML;
                            if ss == 1
                                avg_times = avg_times + (1/num_R)*out_struct.(['ANAM_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.times;
                            end
                        end
                    end
                    semilogy((1/num_A)*avg_ML)
                    if ss == 1
                        total_run_time = [total_run_time; (1/num_A)*sum(avg_times)];
                        run_time_per_iter = [run_time_per_iter; (1/num_A)*mean(avg_times)];
                    end
                end

            case 'FPI'
                avg_ML = zeros(max_iter + 1, 1);
                if ss == 1
                    avg_times = zeros(max_iter, 1);
                end
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_ML = avg_ML + (1/num_R)*out_struct.FPI.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.fun_ML;
                        if ss == 1
                            avg_times = avg_times + (1/num_R)*out_struct.FPI.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.times;
                        end
                    end
                end
                semilogy((1/num_A)*avg_ML)
                if ss == 1
                    total_run_time = [total_run_time; (1/num_A)*sum(avg_times)];
                    run_time_per_iter = [run_time_per_iter; (1/num_A)*mean(avg_times)];
                end

            case 'PAMP'
                for j = 1:len_s
                    avg_ML = zeros(max_iter + 1, 1);
                    if ss == 1
                        avg_times = zeros(max_iter, 1);
                    end
                    for arr = 1:num_A
                        for r = 1:num_R
                            avg_ML = avg_ML + (1/num_R)*out_struct.(['PAMP_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.fun_ML;
                            if ss == 1
                                avg_times = avg_times + (1/num_R)*out_struct.(['PAMP_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.times;
                            end
                        end
                    end
                    semilogy((1/num_A)*avg_ML)
                    if ss == 1
                        total_run_time = [total_run_time; (1/num_A)*sum(avg_times)];
                        run_time_per_iter = [run_time_per_iter; (1/num_A)*mean(avg_times)];
                    end
                end

            case 'SDP'
                avg_ML = zeros(max_iter + 1, 1);
                if ss == 1
                    avg_times = 0;
                end
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_ML = avg_ML + (1/num_R)*out_struct.SDP.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.fun_ML;
                        if ss == 1
                            avg_times = avg_times + (1/num_R)*out_struct.SDP.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.times;
                        end
                    end
                end
                semilogy((1/num_A)*avg_ML)
                if ss == 1
                    total_run_time = [total_run_time; (1/num_A)*avg_times];
                    run_time_per_iter = [run_time_per_iter; Inf];
                end

            case 'WLS'
                avg_ML = zeros(max_iter + 1, 1);
                if ss == 1
                    avg_times = 0;
                end
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_ML = avg_ML + (1/num_R)*out_struct.WLS.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.fun_ML;
                        if ss == 1
                            avg_times = avg_times + (1/num_R)*out_struct.WLS.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.times;
                        end
                    end
                end
                semilogy((1/num_A)*avg_ML)
                if ss == 1
                    total_run_time = [total_run_time; (1/num_A)*avg_times];
                    run_time_per_iter = [run_time_per_iter; Inf];
                end
        end
    end
    set(gca, 'YScale', 'log')
    legend(names_to_legend, 'Interpreter', 'latex', 'FontSize', 14)
    xlim([0 max_iter + 1])
    xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$\textrm{ML\ Function\ Value}$', 'Interpreter', 'latex', 'FontSize', 14)
    title(['$\textrm{Average\ ML\ Function\ Value\ for\ }\sigma=', num2str(sig), '$'], 'Interpreter', 'latex', 'FontSize', 14)
    savefig("output/plots/ML_fun_"+formattedDateTime+"_sigma"+num2str(ss)+"_N"+num2str(N)+"_.fig")
    hold off
end

time_table = table(names_to_legend, total_run_time, run_time_per_iter);
out_struct.time_table = time_table;