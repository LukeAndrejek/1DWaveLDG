function [x,w] = gauss_points_and_weights(n)

%-----------------------------------------------------------------------
%
%  [x,w] = gauss_points_and_weights(n)
%   
%  Returns the Gauss points, x (in ascending order) and their 
%  respective weights, w for an n point Gauss-Legendre quadrature.  
%
%  Written by Ethan Kubatko
%
%-----------------------------------------------------------------------
%
%  Input:
%  ------
%    n:  number of Gauss points
%
%
%  Output:
%  -------
%    x:  1 by n row vector of the n Gauss points
%    w:  1 by n row vector of the respective weights
%
%-----------------------------------------------------------------------
%  NOTE:  An n point Gauss-Legendre quadrature will integrate a 2*n-1
%         degree polynomial exactly.
%-----------------------------------------------------------------------

% Compute the first n+2 legendre polynomials

p = n+1;
P = zeros([p+1,p+1]);
for l = 0:p
    for k = 0:floor(l/2) 
        P(l+1,l-2*k+1) = P(l+1,l-2*k+1) + 1/(2^l)*(-1)^k*factorial(2*l-2*k)...
                         /(factorial(k)*factorial(l-k)*factorial(l-2*k));
    end
end
P = flipdim(P,2);

% Find the n Gauss points

if n==1
    X = 0;
else
    X = roots(P(n+1,:));
end

% Sort the results in ascending order

x = sort(X');

% Compute the corresponding n Gauss weight

w = 2*(1-x.^2)./((n+1)*polyval(P(n+2,:),x)).^2;



