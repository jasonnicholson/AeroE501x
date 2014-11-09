function x = besselzero(n, k, kind)
% besselzero calculates the zeros of Bessel function of the first and second kind 
%
%   x = besselzero(n)
%   x = besselzero(n, k)
%   x = besselzero(n, k, kind)
%
%% Inputs
% * n - The order of the bessel function. n can be a scalar, vector, or
%       matrix.  n can be positive, negative, fractional, or any
%       combinaiton.
% * k - The number of postive zeros to calculate.  When k is not supplied,
%       k = 5 is the default. k must be a scalar.
% * kind - kind is either 1 or 2. When kind is not supplied, default is
%          kind = 1.
%
%% Outputs
% * x - The calculated zeros.  size(x) = [size(n) k].
% 
%% Description
% besselzero calculates the first k positive zeros of nth order bessel
% function of the first or second kind.  Note, that zero is not included as
% the first zero.
%
%% Algorithm
% the first three roots of any order bessel can be approximated by a simple
% equations.  These equations were generated using a least squares fit of
% the roots from orders of n=0:10000. The approximation is used to start
% the iteration of Halley's method.  The 4th and higher roots can be
% approximated by understanding the roots are regularly spaced for a given
% order.  Once the 2nd and 3rd roots are found, the spacing can be
% approximated by the distance between the 2nd and 3rd root.  Then again
% Halley's method can be applied to precisely locate the root.
%
%% Example
%   n = 1;
%   k = 10;
%   kind = 1;
%   z = besselzero(n, k, kind);
%   x = linspace(0, z(end), 1000);
%   y = besselj(n, x);
%   plot(x, y, z, zeros(size(z)),'x')
%

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

%% Input checking
assert(nargin >=1 | nargin <=3,'Wrong number of input arguments.');

% Take care of default cases of k and kind
if nargin < 2
    k = 5;
end

if nargin < 3
    kind = 1;
end

assert(isscalar(kind) & any(kind == [1 2]), '''kind''must be a scalar with value 1 or 2 only.');
assert(isscalar(k) & fix(k)==k & k>0, 'k must a positive scalar integer.');
assert(all(isreal(n(:))), 'n must be a real number.');

%% Setup Arrays

% output size
nSize = size(n);
if nSize(end) ==1
    outputSize = [nSize(1:end-1) k];
else
    outputSize = [nSize k];
end
% number of orders for each kth root
nOrdersPerRoot = prod(outputSize(1:end-1));
x = nan(outputSize);

% negative orders have the same roots as the positive orders
n = abs(n);

%% Solve for Roots
switch kind
    case 1
        % coefficients and exponent are from least squares fitting the k=1,
        % n=0:10000.
        coefficients1j = [2.41383745006265;4.60333884879062;-0.00879177153608845];
        exponent1j = -1.97630722000738;
        % guess for k = 1
        x(((1:nOrdersPerRoot)')') = coefficients1j(1) + coefficients1j(2)*n(:) + coefficients1j(3)*(n(:)+1).^(exponent1j);
        % find first zero
        x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));
        
        if k >= 2
            % coefficients and exponent are from least squares fitting the k=2,
            % n=0:10000.
            coefficients2j = [5.63229987996047;4.60333884784611;-0.116235793367165];
            exponent2j = -0.955326192274869;
            % guess for k = 2
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2j(1) + coefficients2j(2)*n(:) + coefficients2j(3)*(n(:)+1).^(exponent2j);
            % find second zero
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
        end
        
        if k >= 3
            % coefficients and exponent are from least squares fitting the k=3,
            % n=0:10000.
            coefficients3j = [8.85083146233771;4.60333884233918;-0.223115212659848];
            exponent3j = -0.872441208450964;
            % guess for k = 3
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3j(1) + coefficients3j(2)*n(:) + coefficients3j(3)*(n(:)+1).^(exponent3j);
            % find second zero
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
        end
    case 2
        % coefficients and exponent are from least squares fitting the k=1,
        % n=0:10000.
        coefficients1y = [0.804622204593045;4.60333884793812;0.0810993745488260];
        exponent1y = -1.31499999976479;
        % guess for k = 1
        x((1:nOrdersPerRoot)') = coefficients1y(1) + coefficients1y(2)*n(:) + coefficients1y(3)*(n(:)+1).^(exponent1y);
        % find first zero
        x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));
        
        if k >= 2
            % coefficients and exponent are from least squares fitting the k=2,
            % n=0:10000.
            coefficients2y = [4.02305932678643;4.60333884907597;-0.0641172320834957];
            exponent2y = -1.03379588489646;
            % guess for k = 2
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2y(1) + coefficients2y(2)*n(:) + coefficients2y(3)*(n(:)+1).^(exponent2y);
            % find second zero
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
        end
        
        if k >= 3
            % coefficients and exponent are from least squares fitting the k=3,
            % n=0:10000.
            coefficients3y = [7.24155745194042;4.60333884552478;-0.168875347977873];
            exponent3y = -0.907056015534608;
            % guess for k = 3
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3y(1) + coefficients3y(2)*n(:) + coefficients3y(3)*(n(:)+1).^(exponent3y);
            % find second zero
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
        end
    otherwise
        error('Code should never get here.');
end


if k >= 4
    % difference between 2nd and 3rd root
    rootSpacing = x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') - x((nOrdersPerRoot+1:2*nOrdersPerRoot)');
    
    for iRoot = 4:k
        % guesses for remaining roots x[k] = rootSpacing + x[k-1]
        x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = rootSpacing(:) + x(((iRoot-2)*nOrdersPerRoot+1:(iRoot-1)*nOrdersPerRoot)');
        % find the remaining zeros
        x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, k, x0, kind), n, x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)'));
    end
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
end

end