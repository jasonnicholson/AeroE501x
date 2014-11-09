function [cost, residuals, c] = costFunction(exponent, k, n, z)

A = [ones(numel(k),1) n(:) (n(:) + 1).^exponent(1) (n(:) + 2).^exponent(2)];
c = A\z(:);
residuals = reshape(z(:) - A*c, size(z)); 
cost = norm(residuals(:), 2);
