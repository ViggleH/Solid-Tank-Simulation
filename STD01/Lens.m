function [x1,y1,x2,y2,discrim] = Lens(xLine,yLine,m,xEll,yEll,ecc,bEll)

g = yLine - m*(xLine); %g is yintercept from line equation: y = mx + b

aEll = (bEll^2/(1 + ecc^2))^0.5; %a^2 = b^2/(1+ecc^2)

%quadratic coefficients from solving where the equation of a line and
%ellipse are equal

%solving
%(x-h)^2/a^2 + (mx+g-k)^2/b^2 = 1
%bEll = b, aEll = a
%xEll = h, yEll = k

a = (bEll^2 + m^2*aEll^2);

b = 2*(m*aEll^2*(g - yEll) - xEll*bEll^2);

c = xEll^2*bEll^2 + aEll^2*(g - yEll)^2 - aEll^2*bEll^2;

x1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
x2 = (-b + sqrt(b^2 - 4*a*c))/(2*a);

if(x1>x2)
    left=x2;
    x2=x1;
    x1=left;
end

y1 = m*(x1-xLine)+yLine;
y2 = m*(x2-xLine)+yLine;

discrim = b^2 - 4*a*c;

%finding the normal of the ellipse at point x1,y1




return
end