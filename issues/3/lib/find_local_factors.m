function [Lambda0, Lambda_rotated]= find_local_factors(X, r, Lambda0)
% Function to find the sparsest rotation of the leading Principal Components
% of a T by n matrix X. 
% Under sparsity in the loading matrix this will identify the true loading matrix.
%
% returns two arguments:
%   Lambda0: Principal Component estimate
%   Lambda_rotated: Rotation of loading matrix with smallest l1-norm
% See README.txt for more detail
%
% Reference:
% Freyaldenhoven, Simon. "Identification through sparsity in factor models"
% Working paper, 2019.
%addpath('lib')
addpath('source/software/lib')

[T,n] = size(X);
[~,D, V] = svd(X/sqrt(T));
eig_X=diag(D).^2;   
if nargin == 2
    Lambda0=sqrt(n)*V(:,1:r);
end


% compute the rotated solution with minimal l1-norm
[rmat_min] = find_min_rotation(Lambda0); %Finds solution for each point in grid
[~,~,~,Lambda_rotated]=collate_solutions(rmat_min,Lambda0,eig_X); %Combine into candidates
