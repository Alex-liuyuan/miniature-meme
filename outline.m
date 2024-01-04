function x_baseline_clean = outline(x)
    [row, col] = size(x);
    mean_x = mean(x, 1);
    std_x = std(x, 0, 1);
    
    mask = (x > mean_x + 3 * std_x) | (x < mean_x - 3 * std_x);
    x(mask) = NaN;
    
    nan_cols = any(isnan(x));
    x(:, nan_cols) = fillmissing(x(:, nan_cols), 'linear', 'EndValues', 'nearest');

    x_baseline_clean = x;
end