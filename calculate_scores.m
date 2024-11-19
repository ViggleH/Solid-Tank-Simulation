function [score1, score2, score3] = calculate_scores(geo, xInts, yInts, IntMatrix, rayAng)
    % Final function to calculate all four scores
    score1 = calculate_beam_uniformity(IntMatrix); % Beam Uniformity
    score2 = calculate_magnification(IntMatrix); % Magnification
    score3 = calculate_effective_radius(geo, xInts, yInts); % Effective Radius
end
