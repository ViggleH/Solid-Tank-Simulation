function [score, score_01, score_02, score_03] = Solid_Tank_Problem(x, option)

x = STD_proj(x);
eps1 = 0.5;
eps2 = 0.01;
score = zeros(11, 1);
score_01 = zeros(11, 1);
score_02 = zeros(11, 1);
score_03 = zeros(11, 1);

X = [x(1) x(1)+eps1 x(1)-eps1 x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1);
    x(2) x(2) x(2) x(2)+eps1 x(2)-eps1 x(2) x(2) x(2) x(2) x(2) x(2);
    x(3) x(3) x(3) x(3) x(3) x(3)+eps1 x(3)-eps1 x(3) x(3) x(3) x(3);
    x(4) x(4) x(4) x(4) x(4) x(4) x(4) x(4)+eps1 x(4)-eps1 x(4) x(4);
    x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5)+eps2 x(5)-eps2];

parfor i = 1:11
    [score(i), score_01(i), score_02(i), score_03(i)] = Solid_Tank_Sim([X(1:4,i); X(5, i)], option);
end

[score, I] = min(score);
score = -score;
score_01 = -score_01(I);
score_02 = -score_02(I);
score_03 = -score_03(I);

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
