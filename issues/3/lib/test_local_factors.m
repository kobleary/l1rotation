function [has_local_factors, n_small, gamma_n, h_n, Lambda] = test_local_factors(X,r, Lambda, eig_X, alpha_gamma)
% Function to test whether X has local factors
%
% returns up to four arguments:
%    has_local_factors: Logical equal to one if local factors are present
%    n_small: Number of small loadings in sparse rotation
%    gamma_n: Critical value to compare n_small to.
%    Lambda: Rotation of loading matrix with smallest l1-norm
% See README.txt for more detail
%
% Reference:
% Freyaldenhoven, Simon. "Identification through sparsity in factor models"
% Working paper, 2019.

%% Preliminaries
[T,n] = size(X);
if (nargin==2 || isempty(Lambda))
        [~, Lambda]=find_local_factors(X, r);
end
if isempty(X)
    [n,r]=size(Lambda);
    if (round(diag(Lambda'*Lambda)) ~= ones(r,1)*n)
        error("Loading matrix may not be properly normalized. Consider only passing two arguments (X,r).")
    end
    if isempty(eig_X)
        error("If no data X is supplied, at least eigenvalues needed to determine critical values")
    end
else
     eig_X=sort(eig(X'*X/T),'descend');
end
if (nargin<5 || isempty(alpha_gamma))
        alpha_gamma=0.05;
end

c_gamma=-sqrt(2)*erfcinv(2*(1-alpha_gamma/2)); %Identical to norminv(1-alpha_gamma/2)
gamma0=0.03;

h_n=1/(log(n));
expected_small = 1/2*erfc(-h_n/sqrt(2)) - 1/2*erfc(h_n/sqrt(2));%Uses that normcdf(x) = 1/2*erfc(-x/sqrt(2))
gamma= gamma0 + expected_small + c_gamma * sqrt((expected_small*(1-expected_small))/n);
gamma_n=floor(gamma*n);

n_small=sum(abs(Lambda)<h_n, 1);
most_small=sort(n_small,'descend');
has_local_factors= most_small(1) > gamma_n;

