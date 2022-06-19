function output = mib(data, phi_measure, tau, extrapolate)

% -----------------------------------------------------------------------
% This code performs a brute-force search through all possible bipartitions
% of your data, and finds the bipartition that minimizes normalized
% integrated information. If you set extrapolate to 1, this code can find
% an estimate of the ground-truth MIB of a 16-node system within a few hours.
% If you set extrapolate to 0, this code can run for a system of up to
% around 20 nodes within a reasonable amount of time (but not more than
% 20 nodes)
% -----------------------------------------------------------------------


% Generate all possible bipartitions of the network.
all_bipartitions = all_possible_bipartitions(size(data,1));

% Loop through all bipartitions, computing phi for each,
% and take the smallest result (indicating that that result is across the MIB)
if extrapolate
   for i = 1:size(all_bipartitions,1)
        out = extrapolated_phi_calc(data, all_bipartitions(i,:), phi_measure,tau); 
        min_phi(i) = out.phi;
        min_phi_norm(i) = out.phi_norm;
    end
else
    nobs = size(data,2); % number of observations/time steps
    t_range1 = 1: 1: nobs-tau;
    t_range2 = 1+tau: 1: nobs;
    Cov_X = data(:,t_range2)*data(:,t_range2)'/(nobs-tau-1); % equal-time covariance matrix at present
    Cov_Y = data(:,t_range1)*data(:,t_range1)'/(nobs-tau-1); % equal-time covariance matrix at past
    Cov_XY = (data(:,t_range2)*data(:,t_range1)')/(nobs-tau-1); % cross-covariance between the past and present
    for i = 1:size(all_bipartitions,1)
        out = phi_calc(phi_measure,Cov_X,Cov_Y,Cov_XY,all_bipartitions(i,:));
        min_phi(i) = out.phi;
        min_phi_norm(i) = out.phi_norm;
    end
end
phi_ind = find(min_phi_norm==min(min_phi_norm)); % phi_star across the MIB
output.partition=all_bipartitions(phi_ind,:);
output.phi=min_phi(phi_ind);
output.phi_norm=min_phi_norm(phi_ind);