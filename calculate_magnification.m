function score2 = calculate_magnification(IntMatrix)
    % Compute magnification score
    meanIntensity = mean(IntMatrix(:));
    significantDetectors = sum(IntMatrix(:) > meanIntensity * 0.1); % Count significant rays
    magCountRatio = significantDetectors / numel(IntMatrix);
    
    % Scale the score
    if magCountRatio > 0.325
        score2 = (1 / 0.675) * magCountRatio + (1 - (1 / 0.675));
    else
        score2 = 0;
    end
end
