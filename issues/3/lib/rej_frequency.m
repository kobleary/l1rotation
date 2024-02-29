function [rej_freq] = rej_frequency(DGP, n_zeros, alpha_gamma, no_sim)

[T, n, rho_t, rho_i]=readvars('baseline_dgp_param.csv');
has_local_factors=zeros(no_sim,1);

for sim=1:no_sim
    X = create_data(DGP, round(n_zeros), T, n, rho_t, rho_i);
    [has_local_factors(sim), ~, ~] = test_local_factors(X,2,[], [], alpha_gamma);
end

rej_freq=mean(has_local_factors);
end



function X = create_data(DGP, n_zeros, T, n, rho_t, rho_i)

    if DGP == "2d_local_normal"
        lambda_option=1; 
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
    elseif DGP == "2d_local_uniform"
        lambda_option=2; 
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
    elseif DGP == "2d_global_normal"
        lambda_option=1; 
        X = simulate_data_2d(T,n,n-n_zeros,n,rho_t,rho_i,lambda_option);  
    elseif DGP == "2d_local_normal0"
        lambda_option=3; 
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option);
    elseif DGP == "2d_local_normal150"
        lambda_option=1;
        n=150;
        T=n;
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
     elseif DGP == "2d_local_normal200"
        lambda_option=1;
        n=200;
        T=n;
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
     elseif DGP == "2d_local_normal300"
        lambda_option=1;
        n=300;
        T=n;
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
     elseif DGP == "2d_local_normal500"
        lambda_option=1;
        n=500;
        T=n;
        X = simulate_data_2d(T,n-n_zeros,n-n_zeros,n,rho_t,rho_i,lambda_option); 
    else
        error('invalid DGP')
    end

end
