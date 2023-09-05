sig = sigma(ss); 
name_sig = ['sigma',strrep(num2str(sig), '.', '')];
out_struct.(name_met).(name_arr).(name_sig).realizations = cell(num_R,1);