name_arr = ['array',num2str(arr)];
array_to_alg = struct; 
array_to_alg.array_info = arrays_info.general; 
array_to_alg.sensors = arrays.(name_arr).sensors; 
array_to_alg.p = arrays_info.(name_arr).p; 
array_to_alg.bias = arrays_info.(name_arr).bias;
init_to_alg.s0 = start_point_convex_hull*mean(array_to_alg.sensors, 2) + (1 - start_point_convex_hull)*0.5*ones(n, 1);