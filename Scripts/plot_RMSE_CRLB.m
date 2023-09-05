function plot_RMSE_CRLB(n, sigma, len_met, len_m, len_s, num_A, num_R, s_inner_array, methods_to_run, out_struct, arrays_info, names_to_legend,N, formattedDateTime)

% Calculate empirical mean

avg_output = zeros(n+1, length(sigma), len_met, num_A);

for arr = 1:num_A
    for ss = 1:length(sigma)
        met_count = 0; sig = sigma(ss);
        for m = 1:len_m
            met = methods_to_run{m};
            switch met
                case 'ANAM'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        for r = 1:num_R
                            avg_output(:, ss, met_count, arr) = avg_output(:, ss, met_count, arr) + (1/num_R)*out_struct.(['ANAM_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset;
                        end
                    end
                case 'FPI'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        avg_output(:, ss, met_count, arr) = avg_output(:, ss, met_count, arr) + (1/num_R)*out_struct.FPI.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset;
                    end
                case 'PAMP'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        for r = 1:num_R
                            avg_output(:, ss, met_count, arr) = avg_output(:, ss, met_count, arr) + (1/num_R)*out_struct.(['PAMP_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset;
                        end
                    end
                case 'SDP'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        avg_output(:, ss, met_count, arr) = avg_output(:, ss, met_count, arr) + (1/num_R)*out_struct.SDP.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset;
                    end
                case 'WLS'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        avg_output(:, ss, met_count, arr) = avg_output(:, ss, met_count, arr) + (1/num_R)*out_struct.WLS.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset;
                    end
            end
        end
    end
end


% Calculate matrices for RMSE

avg_matrices = zeros(n+1, n+1, length(sigma), len_met, arr);

for arr = 1:num_A
    for ss = 1:length(sigma)
        met_count = 0; sig = sigma(ss);
        for m = 1:len_m
            met = methods_to_run{m};
            switch met
                case 'ANAM'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        for r = 1:num_R
                            temp_vector = out_struct.(['ANAM_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset - avg_output(:, ss, met_count, arr);
                            avg_matrices(:, :, ss, met_count, arr) = avg_matrices(:, :, ss, met_count, arr) + (1/num_R)*(temp_vector*temp_vector');
                        end
                    end
                case 'FPI'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        temp_vector = out_struct.FPI.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset - avg_output(:, ss, met_count, arr);
                        avg_matrices(:, :, ss, met_count, arr) = avg_matrices(:, :, ss, met_count, arr) + (1/num_R)*(temp_vector*temp_vector');
                    end
                case 'PAMP'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        for r = 1:num_R
                            temp_vector = out_struct.(['PAMP_s', num2str(s_inner_array(j))]).(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset - avg_output(:, ss, met_count, arr);
                            avg_matrices(:, :, ss, met_count, arr) = avg_matrices(:, :, ss, met_count, arr) + (1/num_R)*(temp_vector*temp_vector');
                        end
                    end
                case 'SDP'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        temp_vector = out_struct.SDP.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset - avg_output(:, ss, met_count, arr);
                        avg_matrices(:, :, ss, met_count, arr) = avg_matrices(:, :, ss, met_count, arr) + (1/num_R)*(temp_vector*temp_vector');
                    end
                case 'WLS'
                    met_count = met_count + 1;
                    for r = 1:num_R
                        temp_vector = out_struct.WLS.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).realizations{r}.location_and_offset - avg_output(:, ss, met_count, arr);
                        avg_matrices(:, :, ss, met_count, arr) = avg_matrices(:, :, ss, met_count, arr) + (1/num_R)*(temp_vector*temp_vector');
                    end
            end
        end
    end
end


% Calculate average estimated RMSE

RMSE_line = zeros(len_m, length(sigma), arr);
CRLB_line = zeros(length(sigma), 1);

for arr = 1:num_A
    for ss = 1:length(sigma)
        met_count = 0; sig = sigma(ss);
        CRLB_line(ss) = CRLB_line(ss) + (1/num_A)*(sqrt(1/num_R))*arrays_info.(['array', num2str(arr)]).(['sigma',strrep(num2str(sig), '.', '')]).CRLB;
        for m = 1:len_m
            met = methods_to_run{m};
            switch met
                case 'ANAM'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        RMSE_line(met_count, ss, arr) = sqrt(trace(avg_matrices(:, :, ss, met_count, arr)));
                    end
                case 'FPI'
                    met_count = met_count + 1;
                    RMSE_line(met_count, ss, arr) = sqrt(trace(avg_matrices(:, :, ss, met_count, arr)));
                case 'PAMP'
                    for j = 1:len_s
                        met_count = met_count + 1;
                        RMSE_line(met_count, ss, arr) = sqrt(trace(avg_matrices(:, :, ss, met_count, arr)));
                    end
                case 'SDP'
                    met_count = met_count + 1;
                    RMSE_line(met_count, ss, arr) = sqrt(trace(avg_matrices(:, :, ss, met_count, arr)));
                case 'WLS'
                    met_count = met_count + 1;
                    RMSE_line(met_count, ss, arr) = sqrt(trace(avg_matrices(:, :, ss, met_count, arr)));
            end
        end
    end
end

figure(200); hold on
plot(log10(sigma),log10(mean(RMSE_line, 3)'))
plot(log10(sigma),log10(CRLB_line), 'black')
legend([names_to_legend; 'CRLB'], 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14)
grid on; xlabel('$\log_{10}(\sigma)$', 'Interpreter', 'latex');
title('$\mathrm{Average\ Estimated\ RMSE\ (}\log_{10}\mathrm{\ scale)}$')
savefig("output/plots/RMSE_"+formattedDateTime+"_N"+num2str(N)+"_.fig")
hold off