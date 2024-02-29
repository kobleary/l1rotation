function [f] = objectivefcn_spherical(theta,Lambda)
% Calculate objective f Lagrange
[n,r]=size(Lambda);
R=zeros(r,1);

R(1)=cos(theta(1));
for kk=2:r-1
    R(kk)=prod(sin(theta(1:kk-1)),1)*cos(theta(kk));
end
R(r)=prod(sin(theta)');

f = sum(sum(abs(Lambda*R)));
