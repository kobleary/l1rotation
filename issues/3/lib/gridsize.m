function [no_randomgrid] = gridsize(factorno)
% defines number of random draws to start search for local minma from
if factorno==2
        no_randomgrid=300; 
elseif factorno==3
        no_randomgrid=500; 
elseif factorno==4
        no_randomgrid=1000; 
elseif factorno==5
        no_randomgrid=2000; 
elseif (factorno>5 && factorno<9)
        no_randomgrid=3000;
else
        no_randomgrid=5000;
end
