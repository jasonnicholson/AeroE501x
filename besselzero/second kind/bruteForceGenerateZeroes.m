clc; clear; close all;
k = 200;
n = 10000;
z = nan(k, n+1);

initialGuess = 0.800;
stepSize = pi/2;
currentGuess = initialGuess;
nextGuess = initialGuess +stepSize;
yCurrent = bessely(0,currentGuess);
yNext = bessely(0,nextGuess);

for iN = 0:n
    for iK = 1:k    
        % Search for sign change
        while sign(yCurrent) == sign(yNext)
            currentGuess = nextGuess;
            nextGuess = currentGuess + stepSize;
            yCurrent = yNext;
            yNext = bessely(iN, nextGuess);
        end
        
        % find zero in interval
        zero = fzero(@(x) bessely(iN,x), [currentGuess nextGuess]);
        z(iK, iN+1) = zero;
        
        % reset for next loop
        currentGuess = zero + stepSize;
        nextGuess = currentGuess + stepSize;
        yCurrent = bessely(iN,currentGuess);
        yNext = bessely(iN,nextGuess);
    end
    % reset for next loop
        currentGuess = z(1,iN+1) + stepSize;
        nextGuess = currentGuess + stepSize;
        yCurrent = bessely(iN+1,currentGuess);
        yNext = bessely(iN+1,nextGuess);
        disp(iN);
end
[kk, nn] = ndgrid(1:k, 0:n);
zz = z;
save('zeros.mat', 'kk', 'nn', 'zz');