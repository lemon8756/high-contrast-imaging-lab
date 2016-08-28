% zernike_wyant.m

% Created by:   Aaron James Lemmer
% Created on:   Mar 17 2014

% Creates the Zernike polynomial with Wyant sequential ordering j;
% e.g. j = 0 => piston.
% First Wyant index is j = 0, not 1.
        
function Z = zernike_wyant(j, rho, theta)

n = floor(sqrt(j));
m = round((n.^2 + 2*n -j)/2);

if mod(n.^2 + 2*n - j,2) == 0
    %even zernike polynomial
    Z = zrf(n,m,rho).*cos(m*theta);
else %odd zernike polynomial
    Z = zrf(n,m,rho).*sin(m*theta);
end

return

%Zernike radial function (Wyant)
function R = zrf(n, m, rho)
R = 0;
for s = 0:(n-m)
    num = (-1)^s*gamma(2*n-m-s+1);
    denom = gamma(s+1)*gamma(n-s+1)*gamma(n-m-s+1);
    R = R + num/denom*rho.^(2*(n-s)-m);
end