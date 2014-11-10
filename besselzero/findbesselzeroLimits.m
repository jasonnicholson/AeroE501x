% This script explorers the limits of besselzero

clc;clear;close('all');

%% First Kind
upper = 625000;
n = upper;
lower = 0;

while (upper - lower) > eps(n)
    z = besselzero(n, 10,1);
    x = linspace(2*z(1)-z(2), 2*z(10) - z(9),1000);
    y = besselj(n, x);
    figureHandle = figure;
    plot(x,y,z, besselj(n,z),'rx');
    answer = input('u for upper or l for lower: ','s');
    if strcmpi(answer,'u')
        upper = n;
        n = (upper + lower)/2;
    elseif strcmpi(answer, 'l');
        lower = n;
        n = (upper + lower)/2;
    else
        % Do nothing
    end
    close(figureHandle(ishghandle(figureHandle)));
end
nj = n;
fprintf('max n value for 1st kind Bessel function: %1.17e\n', nj);

%% Second Kind
upper = 625000;
n = upper;
lower = 0;

while (upper - lower) > eps(n)
    z = besselzero(n, 10,2);
    x = linspace(2*z(1)-z(2), 2*z(10) - z(9),1000);
    y = bessely(n, x);
    figureHandle = figure;
    plot(x,y,z, bessely(n,z),'rx');
    answer = input('u for upper or l for lower: ','s');
    if strcmpi(answer,'u')
        upper = n;
        n = (upper + lower)/2;
    elseif strcmpi(answer, 'l');
        lower = n;
        n = (upper + lower)/2;
    else
        % Do nothing
    end
    close(figureHandle(ishghandle(figureHandle)));
end
ny = n;
fprintf('max n value for 2nd kind Bessel function: %1.17e\n', ny);