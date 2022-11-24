%Author: Dominic (Zhongda) Huang
%Date: 2022.08.18
%Input: 8 by 1 vector x
%Output: the projection of x on the constraint of STD problem

function p = STD_proj(x)

%x_bl, x_bc, x_lp, x_ma, x_be, bell2, ecc2, d3
constraint = [0, 400; -40 40; 40 100; 40 160; 0 1; 70 120; 0 2.5; 0 400];

%bore radius 52mm + 5mm safeguard
l = 57;

p = x;

for i = 3:(size(x,1)-1)
    if x(i) < constraint(i, 1)
        p(i) = constraint(i, 1);
    else if x(i) > constraint(i, 2)
            p(i) = constraint(i, 2);
        end
    end
end

if(x(1) > 400)
    p(1) = 400;
end

c1 = (x(1) >= 2*l - 2*x(2));
c2 = (x(1) >= 2*l + 2*x(2));

if(c1)
    if(c2)
        if x(2) < constraint(2, 1)
            p(2) = constraint(2, 1);
        else if x(2) > constraint(2, 2)
                p(2) = constraint(2, 2);
            end
        end
    else
        if(x(1) >= -0.5*x(2) + 2*l + 100)
            p(1) = 2*l + 80;
            p(2) = 40;
        else
            if(x(1) >= -0.5*x(2) + 2*l)
            dummy = (dot([x(2); x(1)-2*l], [40, 80])/ dot([40; 80], [40, 80])) * ([40; 80]);
            dummy(2) = dummy(2) + 2*l;
            p(1) = dummy(2);
            p(2) = dummy(1);
            else
                p(1) = 2*l;
                p(2) = 0;
            end
        end
    end
else
    if(c2)
        if(x(1) >= 0.5*x(2) + 2*l + 100)
            p(1) = 2*l + 80;
            p(2) = -40;
        else
            if(x(1) >= 0.5 * x(2) + 2*l)
            dummy = (dot([x(2); x(1)-2*l], [-40, 80])/ dot([-40; 80], [-40, 80])) * ([-40; 80]);
            dummy(2) = dummy(2) + 2*l;
            p(1) = dummy(2);
            p(2) = dummy(1);
            else
                p(1) = 2*l;
                p(2) = 0;
            end
        end
    else
        p(1) = 2 * l;
        p(2) = 0;
    end
end

if (x(8) < 0)
    p(8) = 0;
else
    if (x(8) > 400 - p(1))
        p(8) = 400 - p(1);
    end
end

end