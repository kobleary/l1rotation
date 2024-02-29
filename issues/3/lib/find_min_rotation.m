function [R,fval,exitflag] = find_min_rotation(Lambda)
% find rotations of Lambda with smallest l1-norm across grid of starting
% points
[~,r]=size(Lambda);
no_draws      = gridsize(r); 

fval=zeros(1,no_draws);
exitflag=zeros(1,no_draws);
    
% Create starting points for algorithm
initial_draws=randn(r, no_draws);
initial_draws=normalize(initial_draws, 1, 'norm', 2);

%Convert to polar coordinates
theta=zeros(r-1, no_draws);
for kk=1:r-2
    theta(kk,:)=acot(initial_draws(kk,:)./vecnorm(initial_draws(kk+1:end,:)));
end
theta(r-1,:)=2*acot( (initial_draws(end-1,:) + vecnorm(initial_draws(end-1:end,:))) ./ initial_draws(end,:)); 

%Optimization in polar cooridnates happens w.r.t. theta
angles=theta;

parfor rep=1:no_draws

        starting_point =theta(:, rep);   
        fun2=@(theta) objectivefcn_spherical(theta,Lambda);
        options = optimset('Display','off'); 

        [angles(:,rep),fval(rep),exitflag(rep)] = fminsearch(fun2,starting_point, options); 
end

%Convert back to cartesian coordinates.
R=zeros(r, no_draws);
R(1,:)=cos(angles(1,:));
for kk=2:r-1
    R(kk,:)=prod(sin(angles(1:kk-1,:)),1).*cos(angles(kk,:));
end
R(r,:)=prod(sin(angles),1);   

end

