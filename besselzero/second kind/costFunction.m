function [cost, residuals, c] = costFunction(exponent, k, n, z)

A = [ones(numel(k),1) n(:) (n(:)+1).^exponent];

c = A\z(:);
residuals = z(:) - A*c;
cost = norm(residuals,1);
