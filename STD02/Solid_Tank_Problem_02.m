function [score, score_01, score_02, score_03, score_04] = Solid_Tank_Problem_02(x, gel_option)

x = STD_proj(x);
eps1 = 0.5;
eps2 = 0.01;
score = zeros(17, 1);
score_01 = zeros(17, 1);
score_02 = zeros(17, 1);
score_03 = zeros(17, 1);
score_04 = zeros(17, 1);

X = [x(1) x(1)+eps1 x(1)-eps1 x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1) x(1);
    x(2) x(2) x(2) x(2)+eps1 x(2)-eps1 x(2) x(2) x(2) x(2) x(2) x(2) x(2) x(2) x(2) x(2) x(2) x(2);
    x(3) x(3) x(3) x(3) x(3) x(3)+eps1 x(3)-eps1 x(3) x(3) x(3) x(3) x(3) x(3) x(3) x(3) x(3) x(3);
    x(4) x(4) x(4) x(4) x(4) x(4) x(4) x(4)+eps1 x(4)-eps1 x(4) x(4) x(4) x(4) x(4) x(4) x(4) x(4);
    x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5) x(5)+eps2 x(5)-eps2 x(5) x(5) x(5) x(5) x(5) x(5);
    x(6) x(6) x(6) x(6) x(6) x(6) x(6) x(6) x(6) x(6) x(6) x(6)+eps1 x(6)-eps1 x(6) x(6) x(6) x(6);
    x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7) x(7)+eps2 x(7)-eps2 x(7) x(7);
    x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8) x(8)+eps1 x(8)-eps1];
x01 = X(1, :);
x02 = X(2, :);
x03 = X(3, :);
x04 = X(4, :);
x05 = X(5, :);
x06 = X(6, :);
x07 = X(7, :);
x08 = X(8, :);

parfor i = 1:17
    x_0 = [x01(i); x02(i); x03(i); x04(i); x05(i); x06(i); x07(i); x08(i)]
    [score(i), score_01(i), score_02(i), score_03(i), score_04(i)] = Solid_Tank_Sim_Andy2022(x_0, gel_option);
end

[score, I] = min(score);
score = -score;
score_01 = -score_01(I);
score_02 = -score_02(I);
score_03 = -score_03(I);
score_04 = -score_04(I);

end
