% L = x(1);
% h = x(2);
% dlaser = x(3);
% bEll = x(4);
% ecc = x(5);
% bEll2 = x(6);
% ecc2 = x(7);
% d3 = x(8);

x = [290, 40, 41, 61, 0, 65, 1, 200];

[a,b,c,d,e] = Solid_Tank_Sim_Andy2022(x)

%should return totalscore a = -0.0904

%mine takes about 0.9s to do each run, so if you're using 4 or more cores
%you can expect ~100,000 trials in about 6.25 hours or less.