% Circ.m

% Created by: Aaron James Lemmer
% Created on: May 28 2012

% Creates a circular mask with aperture size D, when given 2D matrices of X
% and Y coordinates.

function z = Circ(x, y, D)
    r = sqrt(x.^2 + y.^2);
    z = double(r<D/2);
    z(r==D/2) = 0.5;