clear;

e = 93;

scores = zeros(1, e);
optimum = zeros(8, e);
process_time = 0;

for k  = 1:e
    load("Water_DQL_Result_" + num2str(k) + ".mat");
    scores(k) = Point_List(x_best_index(end)).Value;
    optimum(:, k) = Point_List(x_best_index(end)).Point;
    process_time = process_time + time;
end
%constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];

%{
for i = 1:e
    for j = 1:5
        scaled_optimum(j, i)  = constraint(j, 1) + optimum(j, i) * (constraint(j, 2) - constraint(j, 1));
    end
    scaled_optimum(:, i) = proj(scaled_optimum(:, i), constraint);
    scaled_optimum(5, i) = scaled_optimum(5, i)  ;
end
%}

[M, I] = min(scores);
-M
x = STD_proj(optimum(:, I))
process_time
[score, score_01, score_02, score_03, score_04] = Solid_Tank_Problem_02(optimum(:, I), 3)

%scaled_optimum(:, I)