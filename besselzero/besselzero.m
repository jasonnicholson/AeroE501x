function x = besselzero(n, k, kind)
% besselzero calculates the zeros of Bessel function of the first and second kind
%
%   x = besselzero(n)
%   x = besselzero(n, k)
%   x = besselzero(n, k, kind)
%
%% Inputs
% * n - The order of the bessel function.
% * k - The number of postive zeros to calculate.  When k is not supplied,
%       k = 5 is the default.
% * kind - kind is either 1 or 2. When kind is not supplied, kind = 1.
%
%% Outputs
% * x - The calculated zeros.  x has length of k.
% Find first k positive zeros of the Bessel function J(n,x) or Y(n,x)
% using Halley's method.
%
%% Description
% besselzero calculates the first k positive zeros of nth order bessel
% function of the first or second kind.
%
%% Example
%   n = 1;
%   k = 10;
%   kind = 1;
%   z = besselzero(n, k, kind);
%   x = linspace(0, z(end), 1000);
%   y = besselj(n, x);
%   plot(x, y, z, zeros(size(z)),'x')

% Originally written by
% Written by: Greg von Winckel - 01/25/05
% Contact: gregvw(at)chtm(dot)unm(dot)edu
%
% Modified and Documented by
% Jason Nicholson 2014-Nov-06
% Contact: jashale@yahoo.com

%% Change Log
% * Original release. 2005-Jan-25, Greg von Winckel.
% * Updated Documentation and commented algorithm. Fixed bug in finding the
%   the first zero of the bessel function of the second kind. 2014-Nov-06,
%   Jason Nicholson.
%

x=nan(k,1);

for iRoot=1:k
    % Initial guess of zeros
    switch kind
        case 1
            switch iRoot
                case 1
                    c = [2.41383745006265;4.60333884879062;-0.00879177153608845];
                    e = -1.97630722000738;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                case 2
                    c = [5.63229987996047;4.60333884784611;-0.116235793367165];
                    e = -0.955326192274869;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                case 3
                    c = [8.85083146233771;4.60333884233918;-0.223115212659848];
                    e = -0.872441208450964;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                otherwise
                    x0 = x(iRoot-1) + x(3)-x(2);
            end
        case 2
            switch iRoot
                case 1
                    c = [0.804622204593045;4.60333884793812;0.0810993745488260];
                    e = -1.31499999976479;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                case 2
                    c = [4.02305932678643;4.60333884907597;-0.0641172320834957];
                    e = -1.03379588489646;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                case 3
                    c = [7.24155745194042;4.60333884552478;-0.168875347977873];
                    e = -0.907056015534608;
                    x0 = c(1) + c(2)*n+c(3)*(n+1).^(e);
                otherwise
                    x0 = x(iRoot-1) + x(3)-x(2);
            end
        otherwise
            error('Code should never get here.');
            
    end
    
    % Do Halley's method
    x(iRoot)=findzero(n,iRoot,x0,kind);
    
end

end



function x=findzero(n,k,x0,kind)
% Uses Halley's method to find a zero given the starting point x0
% http://en.wikipedia.org/wiki/Halley's_method

ITERATIONS_MAX = 100;       % Maximum number of iteration
TOLERANCE_RELATIVE = 1e4;   % 16-4 = 11 significant digits

% Setup loop
error = 1;
loopCount = 0;
x = 1; % Initialization value only.  It is is not used.

% Begin loop
while abs(error)>eps(x)*TOLERANCE_RELATIVE && loopCount<ITERATIONS_MAX
    
    switch kind
        case 1
            a = besselj(n,x0);
            b = besselj((n+1),x0);
        case 2
            a = bessely(n,x0);
            b = bessely((n+1),x0);
    end
    
    xSquared = x0*x0;
    
    error = 2*a*x0*(n*a-b*x0)/(2*b*b*xSquared-a*b*x0*(4*n+1)+(n*(n+1)+xSquared)*a*a);
    
    % Prepare for next loop
    x=x0-error;
    x0=x;
    loopCount=loopCount+1;
    
end

% Handle maximum iterations
if loopCount>ITERATIONS_MAX-1
    warning('Failed to converge to within relative tolerance of %e for n=%f and k=%d in %d iterations', eps(x)*TOLERANCE_RELATIVE, n, k, ITERATIONS_MAX);
    x=NaN;
end

end