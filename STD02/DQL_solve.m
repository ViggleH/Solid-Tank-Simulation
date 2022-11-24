clear;

 for k = 1:3
    j = 1;
    step = 10;
    eps1 = 0.5;
    eps2 = 0.001;
    eps3 = 10^(-3)/0.5;
    max_calls = 8000;
    constraint = [0, 400; -40 40; 40 100; 40 160; 0 1; 70 120; 0 2.5; 0 400];
    gel = k;
        while max_calls >= 17
            x_0 = zeros(5, 1);
            for i = 2:7
                x_0(i) = constraint(i, 1) + rand * (constraint(i, 2) - constraint(i, 1));
            end
            x_0(1) = 2 * 57 + 2 * abs(x_0(2)) + rand * (constraint(1, 2) - (2 * 57 + 2 * abs(x_0(2))));
            x_0(8) = rand * (400 - x_0(1));
            
            f = @(x) Solid_Tank_Problem_02(x, gel);
            [Point_List, x_best_index, stop, time, s, calls, iter] = DQL(f, x_0, step, eps1, eps2, eps3, max_calls,1,2,3, 1);
            switch gel
                case 1
                    parsave("Water_DQL_Result_" + num2str(j) + ".mat", Point_List, x_best_index, stop, time, s, calls, iter);
                case 2
                    parsave("FlexDos_DQL_Result_" + num2str(j) + ".mat", Point_List, x_best_index, stop, time, s, calls, iter);
                case 3
                    parsave("ClearView_DQL_Result_" + num2str(j) + ".mat", Point_List, x_best_index, stop, time, s, calls, iter);
            end
            max_calls = max_calls - calls;
            j = j + 1;          
        end
 end
function parsave(fname, Point_List, x_best_index, stop, time, s, calls, iter)
  save(fname, 'Point_List', 'x_best_index', 'stop', 'time', 's', 'calls', 'iter');
end
