clear;

rng(0);

step = 10;
eps1 = 0.5;
eps2 = 0.001;
eps3 = 10^(-3)/0.5;
max_calls = 5000;
constraint = [0, 400; -40 40; 40 100; 40 160; 0 1];
gel_option = 3;
j = 1;

while max_calls >= 11
    x_0 = zeros(5, 1);
    for i = 2:5
        x_0(i) = constraint(i, 1) + rand * (constraint(i, 2) - constraint(i, 1));
    end
    x_0(1) = 2 * 57 + 2 * abs(x_0(2)) + rand * (constraint(1, 2) - (2 * 57 + 2 * abs(x_0(2))));
    
    f = @(x) Solid_Tank_Problem(x, gel_option);
    [Point_List, x_best_index, stop, time, s, calls, iter] = DQL(f, x_0, step, eps1, eps2, eps3, max_calls, 1,2,3, 1);
    switch gel_option
        case 1
            save("Water_DQL_Result_" + num2str(j) + ".mat");
        case 2
            save("FlexDos_DQL_Result_" + num2str(j) + ".mat");
        case 3
            save("ClearView_DQL_Result_" + num2str(j) + ".mat");
    end
    max_calls = max_calls - calls;
    j = j + 1;
    
end