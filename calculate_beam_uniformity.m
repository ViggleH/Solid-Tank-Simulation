function score1 = calculate_beam_uniformity(IntMatrix)
    % Compute beam uniformity score
    noZpro = nonzeros(IntMatrix(:)); % Ignore zero intensities
    less60 = sum(noZpro < max(noZpro) * 0.6); % Count rays below 60% of max intensity
    beamUni = less60 / length(noZpro); % Ratio of low-intensity rays
    
    % Calculate the beam uniformity score
    score1 = -0.5 * tanh((beamUni - 0.5) * 2 * pi) + 0.5;
end
