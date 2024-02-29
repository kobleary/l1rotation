function [f] = objectivefcn_spherical_2d(theta,Lambda)
% Calculate objective f Lagrange

f = sum(sum(abs(Lambda*[cos(theta); sin(theta)])));
