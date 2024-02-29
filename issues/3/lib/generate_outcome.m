function [y] = generate_outcome(F, noise, beta, beta2, beta3)

[T,~]=size(F);
    switch nargin
        case 3
            y = F*beta + noise*randn(T,1);
        case 4
            y = F*beta + F.^2*beta2 + noise*randn(T,1);
        case 5
            y = F*beta + F.^2*beta2 + F.^3*beta3 + noise*randn(T,1);
    end
    
end

