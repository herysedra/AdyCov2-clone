%Doubling time analysis

% Read - Jewell

sigma = 1/4;
beta = 1.94;
gamma = 1/1.61;

r = 0.5*(sqrt((sigma + gamma )^2 + 4*sigma*( beta - gamma) ) - (sigma + gamma)   )

t_double = log(2)/r

%%

%Li et al

sigma = 1/5.2;
gamma = 1/2.5;

beta = 2.2*gamma;

r = 0.5*(sqrt((sigma + gamma )^2 + 4*sigma*( beta - gamma) ) - (sigma + gamma)   )

t_double = log(2)/r


