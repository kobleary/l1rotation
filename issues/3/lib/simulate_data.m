function [ X, Lambda, F ] = simulate_data( T, n, rho_F, rho_t, rho_i, signal_to_noise,factorno,option)
    %First create Data with an aproximate factor structure.
    %Set parameters
    
    % option-1 standard
    % option-2 no exact sparsity, but rather an approximate version thereof
    
    %Create common component
    
    if factorno==2
       pk = [round(n^(.9)) round(n^(.85))];
    elseif factorno==3
       pk = [round(n^(.9)) round(n^(.85)) round(n^(.8))];
    elseif factorno==4
       pk = [n round(n^(.9)) round(n^(.8)) round(n^(.75))];
    elseif factorno==5
       pk = [n round(n^(.9)) round(n^(.85)) round(n^(.8)) round(n^(.75))];
    elseif factorno==8
       pk = [n n round(n^(.9)) round(n^(.85)) round(n^(.8)) round(n^(.75)) round(n^(.7)) round(n^(.6))];
    else
        print("invalid number of factors")
    end

    %pk = [n n n n n n n n];
    r = length(pk);
    pk = round(pk);

    Lambda=zeros(n,r);
    index=zeros(n,r);
    for k=1:r
        index(:,k) = randperm(n);
        for j = 1: pk(k)
            Lambda(index(j,k),k) = 1;
        end
    end
    Lambda = Lambda .* (ones(n,r) + randn(n,r));
    
    % add non exact sparsity entries
    if option~=1
        for k=1:r
            for j = (pk(k)+1):n
                Lambda(index(j,k),k) = sqrt(1/n)*randn(1,1);
            end
        end
    end

    %Generate Data
    %F
    MU_F=zeros(T,factorno);
    SIGMA_F=[1.0   rho_F    0.1     0.3     0.1     0.3     0.1     0.3;
             rho_F 1.0      rho_F   0.0     0.0     0.0     0.0     0.0;
             0.1   rho_F    1.0     rho_F   0.0     0.0     0.0     0.0;
             0.3   0.0      rho_F   1.0     rho_F   0.0     0.0     0.0;
             0.1   0.0      0.0     rho_F   1.0     rho_F   0.0     0.0;
             0.3   0.0      0.0     0.0     rho_F   1.0     rho_F   0.0;
             0.1   0.0      0.0     0.0     0.0     rho_F   1.0     rho_F;
             0.3   0.0      0.0     0.0     0.0     0.0     rho_F   1.0;];
    
    SIGMA_F=SIGMA_F(1:factorno,1:factorno);

    F=mvnrnd(MU_F,SIGMA_F);

    %Create idiosyncratic noise

    epsilon = randn(T,n);
    cross_dep = [];
    time_dep=[];
    for i=1:n
        cross_dep(i) = rho_i^(i-1);
    end
    for t=1:T
        time_dep(t) = rho_t^(t-1);
    end

    A = toeplitz(time_dep);
    B = toeplitz(cross_dep);
    e = sqrt(A)* epsilon* sqrt(B);

    %Put stuff together
    X= F * Lambda' + sqrt(signal_to_noise)*e;

end
