clc; clear; close all;
load('zeroes.mat');

optionsMinUnc = optimoptions('fminunc', 'TolX', 5e-9, 'Algorithm', 'quasi-newton','Display','Off');
optionsMinSearch = optimset('TolX', 5e-9, 'Display', 'off');
de = 1.5;

%% k = 1
k = 1;
[e1, fval11, exitflag11, output11] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), [0.3327 0.3171], optionsMinUnc);
disp(output11.message);
[cost1, residuals1, a1] = costFunction(e1, kk(k,:), nn(k,:), zz(k,:));
[e1, fval12, exitflag12, output12] = fminsearch(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), e1, optionsMinSearch);
disp(output12.message);
[exponents11, exponents12] = ndgrid((e1(1)-de):0.1:(e1(1)+de),(e1(2)-de):0.1:(e1(2)+de));
costs1 = arrayfun(@(exponent1, exponent2) costFunction([exponent1 exponent2], kk(k,:), nn(k,:), zz(k,:)), exponents11, exponents12);
figure;
surf(exponents11, exponents12, costs1);
hold('all');
plot3(e1(1),e1(2), cost1, 'xr', 'MarkerSize', 10)
title('k=1');
figure;
semilogy(nn(k,:), abs(residuals1))

%% k = 2
k = 2;
[e2,fval21,exitflag21,output21] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), [0.434 0.3598], optionsMinUnc);
disp(output21.message);
[e2, fval22, exitflag122, output22] = fminsearch(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), e2, optionsMinSearch);
disp(output22.message);
[cost2, residuals2, a2] = costFunction(e2, kk(k,:), nn(k,:), zz(k,:));
[exponents21, exponents22] = ndgrid((e2(1)-de):0.1:(e2(1)+de),(e2(2)-de):0.1:(e2(2)+de));
costs2 = arrayfun(@(exponent1, exponent2) costFunction([exponent1 exponent2], kk(k,:), nn(k,:), zz(k,:)), exponents21, exponents22);
figure;
surf(exponents21, exponents22, costs2);
hold('all');
plot3(e2(1),e2(2), cost1, 'xr', 'MarkerSize', 10)
title('k=2');
figure;
semilogy(nn(k,:), abs(residuals2))

%% k = 3
k = 3;
[e3,fval31,exitflag31,output31] = fminunc(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), [0.434 0.4592], optionsMinUnc);
disp(output31.message);
[e3, fval32, exitflag132, output32] = fminsearch(@(exponent) costFunction(exponent, kk(k,:), nn(k,:), zz(k,:)), e3, optionsMinSearch);
disp(output32.message);
[cost3, residuals3, a3] = costFunction(e3, kk(k,:), nn(k,:), zz(k,:));
[exponents31, exponents32] = ndgrid((e3(1)-de):0.1:(e3(1)+de),(e3(2)-de):0.1:(e3(2)+de));
costs3 = arrayfun(@(exponent1, exponent2) costFunction([exponent1 exponent2], kk(k,:), nn(k,:), zz(k,:)), exponents31, exponents32);
figure;
surf(exponents31, exponents32, costs3);
hold('all');
plot3(e3(1),e3(2), cost3, 'xr', 'MarkerSize', 10)
title('k=3');
figure;
semilogy(nn(k,:), abs(residuals3))