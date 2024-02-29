function [Lambda_0,F_0] = pcaest(X,factorno)
%  Given data X, generate the Principal component estimates \Lambda^0, F^0
%  which only take the first factorno number of factors
    [T,n]=size(X);
    Lambda_0=zeros(n,factorno);          
    %[coeff,score]=pca(X); % compute the PCA loading matrix of the data
    [Fhat,d, V] = svd(X/sqrt(T));
    for k=1:factorno
        Lambda_0(:,k) =  sqrt(n)*V(:,k);
    end
    
    F_0=X*Lambda_0/n;
end

