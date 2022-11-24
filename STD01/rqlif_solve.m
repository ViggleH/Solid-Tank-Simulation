clearvars;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
f = @(x) Solid_Tank_Problem(x);
pts = 5;

[Point_List, x_best_index, time, n] = GridSearch(f, pts, constraint, 1);
save Grid_Search_Result.mat;


clear;

step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 600;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
x_0 = zeros(5, 1);
for i = 1:5
    x_0(i)  = constraint(i, 1) + rand() * (constraint(i, 2) - constraint(i, 1));
end

f = @(x) Solid_Tank_Problem(x);
[Point_List, x_best_index, stop, time, s, calls, iter] = rqlif(f, x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1);
save rqlif_result_01.mat;

clear;

step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 600;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
x_0 = zeros(5, 1);
for i = 1:5
    x_0(i)  = constraint(i, 1) + rand() * (constraint(i, 2) - constraint(i, 1));
end

f = @(x) Solid_Tank_Problem(x);
[Point_List, x_best_index, stop, time, s, calls, iter] = rqlif(f, x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1);
save rqlif_result_02.mat;

clear;

step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 600;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
x_0 = zeros(5, 1);
for i = 1:5
    x_0(i)  = constraint(i, 1) + rand() * (constraint(i, 2) - constraint(i, 1));
end

f = @(x) Solid_Tank_Problem(x);
[Point_List, x_best_index, stop, time, s, calls, iter] = rqlif(f, x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1);
save rqlif_result_03.mat;

clear;

step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 600;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
x_0 = zeros(5, 1);
for i = 1:5
    x_0(i)  = constraint(i, 1) + rand() * (constraint(i, 2) - constraint(i, 1));
end

f = @(x) Solid_Tank_Problem(x);
[Point_List, x_best_index, stop, time, s, calls, iter] = rqlif(f, x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1);
save rqlif_result_04.mat;

clear;

step = 10;
eps1 = 10^(-3);
eps2 = 10^(-12);
eps3 = 10^(-6);
max_calls = 600;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
x_0 = zeros(5, 1);
for i = 1:5
    x_0(i)  = constraint(i, 1) + rand() * (constraint(i, 2) - constraint(i, 1));
end

f = @(x) Solid_Tank_Problem(x);
[Point_List, x_best_index, stop, time, s, calls, iter] = rqlif(f, x_0, step, eps1, eps2, eps3, max_calls, 2, 3, 5, 1);
save rqlif_result_05.mat;