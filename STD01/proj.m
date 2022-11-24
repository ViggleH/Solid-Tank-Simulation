%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: vector x, box-constraints S
%Output: the projection of x on S

function p = proj(x, S)

    p = x;

    for i = 1:size(x,1)
        if x(i) < S(i, 1)
            p(i) = S(i, 1);
        else if x(i) > S(i, 2)
            p(i) = S(i, 2);
            end
        end
    end

end