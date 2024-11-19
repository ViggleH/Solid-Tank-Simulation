function visualize_geometry(geo, xInts, yInts, scale)
    % Function to visualize the tank setup and ray paths with all features
    theta = 0:0.001:2*pi; % Angle range for plotting circles and ellipses

    % Create figure
    figure;
    hold on;
    axis square;
    
    % Plot tank elements: ellipses and circles
    plot_ellipse_and_circles(geo, scale, theta);

    % Plot detector boundaries and central scatter points
    plot_detectors_and_central_points(geo, scale, xInts);

    % Draw ray paths
    plot_ray_paths(xInts, yInts, geo, scale);

    % Finalize visualization
    xlim([-165, 165]);
    ylim([-90, 90]);
    daspect([1, 1, 1]);
    hold off;
end
