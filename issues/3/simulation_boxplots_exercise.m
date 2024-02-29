addpath('lib')

% Set parameter values
rng(1909);
factorno=4;
T=500;
n=300;
rho_F = 0; %AR coefficients in Factors. 
rho_t= 0.3; %Time dependance
rho_i= 0.1; %cross-sectional  dependance
signal_to_noise= 1;
no_sim=500;

tic
%Repeat with exact sparstity and approximate sparsity
for exact_sparsity=1:2
    % Generate matrices to store results
    true_Lambda_normal_all=zeros(n,factorno,no_sim); 
    X_all=zeros(T,n,no_sim); % all simulated data
    F_all=zeros(T,factorno,no_sim); % all Factors
    rmat_min_all=zeros(factorno,[],no_sim); % all rotation column vectors which lead to local minimum of the objective function
    l1_min_all=zeros(no_sim,[]); % local minmum of the objective function
    exitflag_all=zeros(no_sim,[]);
    Lambda_0_all=zeros(n,factorno,no_sim); % pca estimator of loading matrix
    eig_X_all=zeros(n,no_sim);
    
    % Main simulation
    for sim=1:no_sim
        if(exact_sparsity == 1)
            spars = '';
        else
            spars = 'approx_spars_';
        end
        disp(sim)
        % simulate data and find principal component estimates
        [X,true_Lambda,F] = simulate_data( T, n, rho_F, rho_t, rho_i, signal_to_noise, factorno,exact_sparsity); % simulate data
        [Lambda_0,~] = pcaest(X,factorno); % Given data X, generate the principal component estimates \Lambda^0, F^0
        eig_X=sort(eig(X'*X/T),'descend'); % compute the eigenvalues of X'*X
        
        true_Lambda_normal=normalize(true_Lambda, 1, 'norm', 2)*sqrt(n);

        % for each point in large grid, find the local minimum of l1 norma and their corresponding rotation vector rmat_min_col
        [rmat_min,l1_min,exitflag] = find_min_rotation(Lambda_0);
        [rmat_min_unique,l1_min_unique,value_counts,Lambda_rotated]=collate_solutions(rmat_min,Lambda_0, eig_X);

        writematrix(true_Lambda_normal, ['datastore/simulation/true_lambda_', spars, num2str(sim), '.csv'])
        writematrix(Lambda_rotated, ['datastore/simulation/lambda_rotated_', spars, num2str(sim), '.csv'])
        writematrix(X, ['datastore/simulation/X_', spars, num2str(sim), '.csv'])

        %save(['datastore/simulation/boxplot_all_candidates',num2str(sim)],'rmat_min_unique','l1_min_unique','exitflag','true_Lambda_normal','X','Lambda_0','eig_X', 'F');

     end
     
end

toc

%exit;