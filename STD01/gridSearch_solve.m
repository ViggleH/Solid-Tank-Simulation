clearvars;

constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
f = @(x) Solid_Tank_Problem(x);
pts = 5;

[Point_List, x_best_index, time, n] = GridSearch(f, pts, constraint, 1);
save Water_Grid_Search_Result.mat;