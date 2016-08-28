% zernike_noll.m

% Created by:   Aaron James Lemmer
% Created on:   Feb 18 2013
% Source:       Numerical Simulation of Optical Wave Propagation, by Jason
%               D. Schmidt, pp. 68-69.

% Modified on:  Mar 14 2014
	% Fixed a typo.  Added automatic mapping from Noll index j to (m,n).
        
function Z = zernike_noll(j, r, theta)
%Creates the Zernike polynomial with Noll mode index j;
%e.g. j = 1 => piston
%First Noll index is 1, not 0

n = floor(sqrt(2*j-1) + 0.5) - 1;
if mod(n,2) == 0  %n is even
	m = 2*floor((2*j+1-n*(n+1))/4);
else  %n is odd
	m = 2*floor((2*(j+1)-n*(n+1))/4) - 1;
end

if m==0
    Z = sqrt(n+1)*zrf(n,0,r);
else
    if m > 0 %even zernike polynomial
        Z = sqrt(2*(n+1))*zrf(n,m,r).*cos(m*theta);
    else %odd zernike polynomial
        Z = sqrt(2*(n+1))*zrf(n,m,r).*sin(m*theta);
    end
end
return

%Zernike radial function
function R = zrf(n, m, r)
R = 0;
for s = 0:((n-abs(m))/2)
    num = (-1)^s*gamma(n-s+1);
    denom = gamma(s+1)*gamma((n+m)/2-s+1)*gamma((n-m)/2-s+1);
    R = R + num/denom*r.^(n-2*s);
end