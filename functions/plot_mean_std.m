function plot_mean_std(x, meanVals, stdVals, plotColor)
    % Set default color to blue if not provided by the user
    if nargin < 4 || isempty(plotColor)
        plotColor = 'b'; 
    end

    % Ensure column vectors
    x = x(:);
    meanVals = meanVals(:);
    stdVals = stdVals(:);
    
    % Upper and lower bounds
    upper = meanVals + stdVals;
    lower = meanVals - stdVals;
    
    % Create shaded area
    % We apply the chosen color to the fill, and rely on FaceAlpha to make it look lighter
    h = fill([x; flipud(x)], [upper; flipud(lower)], ...
        plotColor, 'EdgeColor', 'none');
    set(h, 'FaceAlpha', 0.3);
    
    hold on;
    
    % Plot mean line
    % We use 'Color' instead of a positional argument to support RGB triplets
    plot(x, meanVals, 'Color', plotColor, 'LineWidth', 2);
    
    hold off;
    
    % Aesthetics
    xlabel('X');
    ylabel('Value');
    title('Mean with Standard Deviation');
    grid off;
end