function [R,fval,sol_frequency,Lambda_rotated]=collate_solutions(rmat_min,Lambda_0,eig_X)
    
    [n,~]=size(Lambda_0);
    [factorno, no_randomgrid]=size(rmat_min); 
    
    epsilon_rot=0.05; %Defines how close columns need to be to be treated as equal
    
    %% Create unique candidates
    l1_min=sum(abs(Lambda_0*rmat_min));
    [l1_min_sort,sort_index] = sort(l1_min); % sort in an ascending order the rotation columns by the l1-norms they generate  
    rmat_min_sort=rmat_min(:,sort_index);
    rmat_min_sort=rmat_min_sort.*repmat(sign(rmat_min_sort(1,:)),factorno,1); % Normalise sign

    % Make "close" columns equal to the ones with the smallest l1-norms in each group.
    for i=1:no_randomgrid
        for j=i:no_randomgrid
            if norm(rmat_min_sort(:,i)-rmat_min_sort(:,j),2)/factorno <epsilon_rot
              rmat_min_sort(:,j)=rmat_min_sort(:,i);
              l1_min_sort(j)=l1_min_sort(i);
            end
        end
    end

    [l1_min_unique,starting_indices,group_id] = unique(l1_min_sort); 
    value_counts = [l1_min_unique', accumarray(group_id,1)]; % count the occurence of each unqiue rotation column vector
    %[value_counts_sorted,index_sort] = sortrows(value_counts,2,'descend'); % Option to sort rotations by frequency
    rmat_min_unique = rmat_min_sort(:,starting_indices);
    
    %% Consolidate candidates
     Lambda_rotated = Lambda_0*rmat_min_sort(:,1); %First candidate
     R=rmat_min_unique(:,1);
     sol_frequency = value_counts(1,2);
     fval =value_counts(1,1);
     
     for kk=2:length(starting_indices)
        temp=[Lambda_rotated Lambda_0*rmat_min_unique(:,kk)]; 
        if (min(eig(temp'*temp)/n)> sqrt(1/factorno)/3 && min(eig(temp'*temp)/n)> min(eig(Lambda_rotated'*Lambda_rotated)/n)/4)
            Lambda_rotated = temp;
            R = [R rmat_min_unique(:,kk)];
            sol_frequency = [sol_frequency value_counts(kk,2)];
            fval =[fval value_counts(kk,1)];
        end
     end
    
    %% If too few candidates, fill in with PCs that are far from collinear
    [~,r_temp]=size(Lambda_rotated);
    I=eye(factorno);
    while r_temp < factorno
        min_eig=zeros(factorno,1);
        for ell=1:factorno
            temp=[Lambda_rotated Lambda_0(:,ell)];
            min_eig(ell)=min(eig(temp'*temp));
        end
        [~,index]=max(min_eig.*sqrt(eig_X(1:factorno))); 
        Lambda_rotated = [Lambda_rotated Lambda_0(:,index)]; 
        R = [R I(:,index)]; 
        fval = [fval sum(abs(Lambda_0(:,index)))]; 
        sol_frequency =[sol_frequency 0]; 
        [~,r_temp]=size(Lambda_rotated);
    end
    