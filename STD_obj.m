function [score1, score2, score3] = STD_obj(x)
    % Main function for calculating and visualizing scores for refractive index bath container
    
    % Constants
    scale = 1;
    numProjections = 720;
    rot = 360;
    laserType = 'fan';
    N = 512;
    
    % Initialize parameters
    rPool = 40; % Ray pool size
    numDet = 2048; % Number of detectors
    
    % Step 1: Setup geometry and refractive parameters
    geo = setup_geometry(rPool, numDet, laserType, N, x);
    
    % Step 2: Compute intersection points and intensities
    [xInts, yInts, IntMatrix, rayAng] = compute_intersections(geo);
    
    % Step 3: Visualize the setup and ray paths
    visualize_geometry(geo, xInts, yInts, scale);
    
    % Step 4: Calculate scores
    [score1, score2, score3] = calculate_scores(geo, xInts, yInts, IntMatrix, rayAng);
end
