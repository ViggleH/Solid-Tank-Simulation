clear;

step = 0.2;
eps1 = 0.0125;
eps2 = 0.0025;
eps3 = 10^(-3)/0.0025;
max_calls = 3125;
constraint = [200, 400; -30 30; 40 100; 40 80; 0 50];
lb = [0; 0; 0; 0; 0];
ub = [1; 1; 1; 1; 1];



f = @(x) Scaled_Solid_Tank_Problem(x);
options = optimoptions('surrogateopt', 'MaxFunctionEvaluations', max_calls, 'MinSampleDistance', eps2);
[x,fval,exitflag,output] = surrogateopt(f,lb,ub,options);
save("Water_SO_Result.mat")
