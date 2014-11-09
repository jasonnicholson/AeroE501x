clc; clear; close all;
load('zeroes.mat');

options = optimoptions('fminunc', 'TolX', 5e-9);

%% k = 1
k = 1;
[e1,fval1,exitflag1,output1] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), -1.96, options);
disp(output1.message);
[cost1, residuals1, a1] = costFunction(e1, kk(k,:), nn(k,:), zz(k,:));
exponents1 = (e1-1):0.05:(e1+1);
costs1 = arrayfun(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), exponents1);
figure;
plot(exponents1, costs1, e1, cost1, 'xr', 'MarkerSize', 10)
title('k=1');
figure;
semilogy(nn(k,:), abs(residuals1))

%% k = 2
k = 2;
[e2,fval2,exitflag2,output2] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), -1.315, options);
disp(output2.message);
[cost2, residuals2, a2] = costFunction(e2, kk(k,:), nn(k,:), zz(k,:));
exponents2 = (e2-1):0.05:(e2+1);
costs2 = arrayfun(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), exponents2);
figure;
plot(exponents2, costs2, e2, cost2, 'xr', 'MarkerSize', 10)
title('k=2');
figure;
semilogy(nn(k,:), abs(residuals2))

%% k = 3
k = 3;
[e3,fval3,exitflag3,output3] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), -1.315, options);
disp(output3.message);
[cost3, residuals3, a3] = costFunction(e3, kk(k,:), nn(k,:), zz(k,:));
exponents3 = (e3-1):0.05:(e3+1);
costs3 = arrayfun(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), exponents3);
figure;
plot(exponents3, costs3, e3, cost3, 'xr', 'MarkerSize', 10)
title('k=3');
figure;
semilogy(nn(k,:), abs(residuals3))