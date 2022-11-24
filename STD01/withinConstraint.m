function a = withinConstraint(x, C)

a = 1;
for i = 1:size(C, 1)
    if (x(i) < (C(i, 1)) || x(i) > (C(i, 2)))
        a = 0;
        break;
    end
end

end