function score3 = calculate_effective_radius(geo, xInts, yInts)
    % Identify the shortest ray that passes through the center of the tank
    shortRay = zeros(size(xInts, 1), 1);
    for i = 1:size(xInts, 1)
        if is_valid_ray(xInts(i, :), yInts(i, :), geo)
            % Compute chord length
            intPts = [xInts(i, 6), yInts(i, 6); xInts(i, 5), yInts(i, 5)];
            shortRay(i) = pdist(intPts, 'euclidean');
        else
            shortRay(i) = NaN;
        end
    end
    
    % Find the effective radius
    shortRay = shortRay(~isnan(shortRay));
    if isempty(shortRay)
        effRad = 0; % No valid rays
    else
        effRad = min(shortRay); % Smallest chord length
    end
    
    % Compute the score
    score3 = 0.5 * tanh((10) * (effRad / geo.r1 - 0.96) * 2 * pi) + 0.5;
end

function valid = is_valid_ray(xInts, yInts, geo)
    % Helper function to check if a ray is valid
    valid = (max(yInts) < (geo.arrayW + geo.k)) && ...
            (min(yInts) > -(geo.arrayW + geo.k)) && ...
            (nnz(xInts) > 9) && ...
            (xInts(1) == geo.det);
end
