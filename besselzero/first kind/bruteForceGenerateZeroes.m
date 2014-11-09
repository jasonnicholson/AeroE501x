clc; clear; close all;
tic;
k = 200;
n = 10000;
z = nan(k, n+1);

initialGuess = 2.2;
stepSize = 0.25;
currentGuess = initialGuess;
nextGuess = initialGuess +0.5;
yCurrent = besselj(0,currentGuess);
yNext = besselj(0,nextGuess);

for iN = 0:n
    for iK = 1:k    
        % Search for sign change
        while sign(yCurrent) == sign(yNext)
            currentGuess = nextGuess;
            nextGuess = currentGuess + stepSize;
            yCurrent = yNext;
            yNext = besselj(iN, nextGuess);
        end
        
        % find zero in interval
        zero = fzero(@(x) besselj(iN,x), [currentGuess nextGuess]);
        z(iK, iN+1) = zero;
        
        % reset for next loop
        currentGuess = zero + stepSize;
        nextGuess = currentGuess + stepSize;
        yCurrent = besselj(iN,currentGuess);
        yNext = besselj(iN,nextGuess);
    end
    % reset for next loop
        currentGuess = z(1,iN+1) + stepSize;
        nextGuess = currentGuess + stepSize;
        yCurrent = besselj(iN+1,currentGuess);
        yNext = besselj(iN+1,nextGuess);
        disp(iN);
end
toc;

[kk, nn] = ndgrid(1:k, 0:n);
zz = z;
save('zeroes.mat', 'kk', 'nn', 'zz');

for iN=0:100:n
    figure('name', sprintf('k=%d',iN),'position', [993,1,927,973]);
    x = linspace(zz(1,iN+1)-2*pi, zz(10,iN+1)+1, 200);
    y = besselj(iN,x);
    plot(x,y,zz(1:10,iN+1),besselj(iN,zz(1:10,iN+1)),'rx','MarkerSize',10);
    axis('tight');
    drawnow();
end