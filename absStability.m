function [meshX, meshY, R_Z] = absStability(A, b, plotSize)
    b = reshape(b, 1, length(b));
    [meshX, meshY] = meshgrid(linspace(-plotSize, plotSize, 20*plotSize), linspace(-plotSize, plotSize, 20*plotSize));
    Z = meshX+meshY*1i;
    I = eye(length(b), length(b));
    e = ones(length(b), 1);
    % TODO: this is badly optimized, but it should be as it needn't be
    % performant
    R = @(z) det(I-z*A+z*e*b)/det(I-z*A);
    R_Z = zeros(size(Z));
    for i = 1:size(Z, 2)^2
        R_Z(i) = R(Z(i));
    end
end