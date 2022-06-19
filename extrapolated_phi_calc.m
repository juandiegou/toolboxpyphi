
function output = extrapolated_phi_calc(data, part, phi_measure,...
    tau)

% -----------------------------------------------------------------------
% This function extrapolates integrated information across a given
% partition to infinite observations. It does so following well-established
% methods for extrapolating entropy and mutual information to infinte
% observations: it computes integrated information as a function of 1/N,
% where N is number of observations. It then fits a line to the resulting
% distribution, and picks the extrapolated information value as the
% y-intercept of this line (i.e. where 1/N=0, which is where N approaches
% infinity)
% -----------------------------------------------------------------------

% Calculate phi using all available observations - this is the value you 
% would get without extrapolation
full_out = phi_comp(data, part, phi_measure,tau,0);

% Store the partition used in the above calculation of phi, so you can use
% the same partition in all the following calculations of phi with fewer
% data points
partition=full_out.partition;
output.partition=full_out.partition;
output.entropy=full_out.entropy;

N=size(data,2);

% Calculate phi for different fractions of all available observations:

% We will split the data into 10 chunks and compute integrated information
% for each of those, then 9 chunks and compute integrated information for
% each of those, etc up to splitting the data in half and computing
% integrated information for each half:
div_vec = 10:-1:2;

for i = 1:length(div_vec)
    div_inds(1:div_vec(i),1) = 1:N/div_vec(i):N;
    div_inds(1:div_vec(i),2) = (1:N/div_vec(i):N)+N/div_vec(i)-1;
    div_inds=round(div_inds);
    
    for d = 1:div_vec(i)
        dat = data(:,div_inds(d,1):div_inds(d,2));
        nobs=size(dat,2);
        timerange1 = 1: nobs-tau;
        timerange2 = 1+tau: nobs;
        Cov_X = dat(:,timerange2)*dat(:,timerange2)'/(nobs-tau-1); % covariance matrix of the present
        Cov_Y = dat(:,timerange1)*dat(:,timerange1)'/(nobs-tau-1); % covariance matrix of the past past
        Cov_XY = (dat(:,timerange2)*dat(:,timerange1)')/(nobs-tau-1); % cross-covariance between the past and present
        out = phi_calc(phi_measure,Cov_X,Cov_Y,Cov_XY,partition);
        phi_norm(d) = out.phi_norm;
        phi(d) = out.phi; 
    end
    phi_norm_y(i)=mean(phi_norm);
    phi_y(i) = mean(phi);
    clear phi phi_norm div_inds 
end

% Fit a line to these sub-sampled estimates:

size_vec = [N./div_vec N]; % N
x=1./size_vec; % 1/N


phi_norm_y(end+1) = full_out.phi_norm; % add what normalized phi was when we took all N observations
phi_norm_y=phi_norm_y';
phi_norm_f = fit(x(~isinf(phi_norm_y) & ~isnan(phi_norm_y))', ...
    phi_norm_y(~isinf(phi_norm_y) & ~isnan(phi_norm_y)), 'poly1'); % fit a line to the observations


% same as above, but for the non-normalized phi values
phi_y(end+1) = full_out.phi;
phi_y = phi_y';
phi_f = fit(x(~isinf(phi_y) & ~isnan(phi_y))', ...
    phi_y(~isinf(phi_y) & ~isnan(phi_y)), 'poly1');


% Store the results
output.phi = phi_f(0); % extrapolated phi value
output.phi_norm = phi_norm_f(0); % extrapolated phi value, normalized
output.x = x; % your 1/n vector
output.phi_y = phi_y; % the phi values corresponding to each point along your x vector
output.phi_f = phi_f; % the line fitted to your values
output.phi_norm_y = phi_norm_y; % the phi (normalized) values corresponding to each point along your x vector
output.phi_norm_f = phi_norm_f; % the line fitted to your normalized values


