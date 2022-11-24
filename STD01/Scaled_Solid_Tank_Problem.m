function score = Scaled_Solid_Tank_Problem(x)
constraint = [200, 400; -30 30; 40 100; 40 80; 0 1];

for i = 1:5
       x(i)  = constraint(i, 1) + x(i) * (constraint(i, 2) - constraint(i, 1));
end

x = proj(x, constraint);
eps1 = 0.5;
eps2 = 0.01;
score = zeros(11, 1);

X = [x(1) x(1)+eps1 x(1)-eps1 x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1);
    x(2) x(2) x(2) x(2)+eps1 x(2)-eps1 x(2) x(2) x(2) x(2) x(2) x(2);
    x(3) x(3) x(3) x(3) x(3) x(3)+eps1 x(3)-eps1 x(3) x(3) x(3) x(3);
    x(4) x(4) x(4) x(4) x(4) x(4) x(4) x(4)+eps1 x(4)-eps1 x(4) x(4);
    x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5)+eps1 x(5)-eps1];

parfor i = 1:11
    [score(i), ~, ~, ~] = Solid_Tank_Sim([X(1:4,i); X(5, i)*1/50]);
end

score = min(score);
score = -score;

%{
if (withinConstraint(x, constraint))
eps = 0.5;
score = zeros(11, 1);

X = [x(1) x(1)+eps x(1)-eps x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1);
    x(2) x(2) x(2) x(2)+eps x(2)-eps x(2) x(2) x(2) x(2) x(2) x(2);
    x(3) x(3) x(3) x(3) x(3) x(3)+eps x(3)-eps x(3) x(3) x(3) x(3);
    x(4) x(4) x(4) x(4) x(4) x(4) x(4) x(4)+eps x(4)-eps x(4) x(4);
    x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5)+eps x(5)-eps];

parfor i = 1:11
    score(i) = Solid_Tank_Sim([X(1:4,i); 1/50 * X(5,i)]);
end

score = min(score);
score = -score;
else
    score = +inf;

end
%}
end