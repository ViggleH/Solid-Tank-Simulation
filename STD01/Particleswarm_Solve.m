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
options = optimoptions('particleswarm', 'MaxFunctionEvaluations', max_calls);
[x,fval,exitflag,output] = particleswarm(f,5,lb,ub,options);
save("Water_ParticalSwarm_Result.mat");
