function plot_Bias(sigma, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, max_iter, names_to_legend,N, formattedDateTime)

for ss = 1:length(sigma)
    sig = sigma(ss);
    figure(100 + ss); hold on
    for m = 1:len_m
        met = methods_to_run{m};
        switch met
            case 'ANAM'
                for j = 1:len_s
                    avg_Bias = zeros(max_iter + 1, 1);
                    for arr = 1:num_A
                        for r = 1:num_R
                            avg_Bias = avg_Bias + (1/num_R)*out_struct.(['ANAM_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.biasV;
                        end
                    end
                    semilogy((1/num_A)*avg_Bias)
                end
            case 'FPI'
                avg_Bias = zeros(max_iter + 1, 1);
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_Bias = avg_Bias + (1/num_R)*out_struct.FPI.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.biasV;
                    end
                end
                semilogy((1/num_A)*avg_Bias)
            case 'PAMP'
                for j = 1:len_s
                    avg_Bias = zeros(max_iter + 1, 1);
                    for arr = 1:num_A
                        for r = 1:num_R
                            avg_Bias = avg_Bias + (1/num_R)*out_struct.(['PAMP_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.biasV;
                        end
                    end
                    semilogy((1/num_A)*avg_Bias)
                end
            case 'SDP'
                avg_Bias = zeros(max_iter + 1, 1);
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_Bias = avg_Bias + (1/num_R)*out_struct.SDP.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.biasV;
                    end
                end
                semilogy((1/num_A)*avg_Bias)
            case 'WLS'
                avg_Bias = zeros(max_iter + 1, 1);
                for arr = 1:num_A
                    for r = 1:num_R
                        avg_Bias = avg_Bias + (1/num_R)*out_struct.WLS.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.biasV;
                    end
                end
                semilogy((1/num_A)*avg_Bias)
        end
    end
    set(gca, 'YScale', 'log')
    legend(names_to_legend, 'Interpreter', 'latex', 'FontSize', 14)
    xlim([0 max_iter + 1])
    xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$\textrm{Bias}$', 'Interpreter', 'latex', 'FontSize', 14)
    title(['$\textrm{Average\ Squared\ Norm\ of\ Bias\ for\ }\sigma=', num2str(sig), '$'], 'Interpreter', 'latex', 'FontSize', 14)
    savefig("output/plots/Bias_"+formattedDateTime+"_sigma"+num2str(ss)+"_N"+num2str(N)+"_.fig")
    hold off
end