function [ X, Lambda, truth_normal, F, r] = simulate_data_2d(T, m1,m2,n,rho_t, rho_i, lambda_option)
    %rho_t, rho_i correlation in error terms.

    % Dep defines the structure of the factors.
    Dep = [1 m1; n+1-m2 n];
    dim=size(Dep);
    r = dim(1);   
    F= randn(T, r);
    F = F*chol([1 0.3; 0.3 1]);% Cholesky decomposition to transform factors to have 0.3 theoretical correlation and variance 1.
    Lambda=zeros(n,r);
    if lambda_option==1
        for k = 1:r
            Lambda(Dep(k,1):Dep(k,2),k) = ones(Dep(k,2)-Dep(k,1)+1,1)+normrnd(0,1,[Dep(k,2)-Dep(k,1)+1,1]); 
        end
    elseif lambda_option==2
        for k = 1:r
            Lambda(Dep(k,1):Dep(k,2),k) = unifrnd(0.1,1.9,[Dep(k,2)-Dep(k,1)+1,1]); 
        end
    elseif lambda_option==3
        for k = 1:r
            Lambda(Dep(k,1):Dep(k,2),k) = normrnd(0,1,[Dep(k,2)-Dep(k,1)+1,1]); 
        end
    end
    
    %Create noise
    epsilon = randn(T,n);
    cross_dep = [];
    time_dep=[];
    for i=1:n
        cross_dep(i) = rho_i^(i-1);
    end
    for t=1:T
        time_dep(t)=rho_t^(t-1);
    end
     
    A = toeplitz(time_dep);
    B = toeplitz(cross_dep);
    e = sqrt(A)* epsilon* sqrt(B);
    %Now create X
    X= F * Lambda' + 1.0*e;
    
    %Also save normalised loading matrix
    truth_normal=[Lambda];
    for k=1:r
        truth_normal(:,k)=[sqrt(n)*Lambda(:,k)/norm(Lambda(:,k))]; % normalize i1 by columns to unit length
    end
end

