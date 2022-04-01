% IDX2RAD Transform TOMS REFL matrix indices to radians.
%
% [theta,phi] = idx2rad(i,j,sy,sx)
%
% $Id: idx2rad.m,v 1.5 2006/05/17 14:39:18 danji Exp $

function [theta,phi] = idx2rad(i,j,sy,sx);
    sy = double(sy);
    sx = double(sx);
    CONST.d2r = pi/180;
    dx = 2*pi/sx;
    dy = pi/sy;
    phi = pi-dy/2-(double(i)-1)*dy;
    theta = (double(j)-1)*dx-pi+dx/2;
return