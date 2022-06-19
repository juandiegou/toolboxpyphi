function output = phi_calc(phi_measure,Cov_X,Cov_Y,Cov_XY, partition)

% -----------------------------------------------------------------------
% This is the mathematical meat of this toolbox. It can calculate either phi_G
% or phi_M (which is the time-reverse of "stochastic interaction") from the
% covariance matrix of the present of your system (Cov_X), the covariance
% matrix of the past of your system (Cov_Y), the cross-covariance between
% the past and present (Cov_XY), and a candidate system partition
% -----------------------------------------------------------------------

output.partition=partition; % store the candidate partition in the output

% set up flags
if strcmp(phi_measure,'phi_M')
    phi_M_calc=1;
    phi_G_calc=0;
elseif strcmp(phi_measure, 'phi_G')
    phi_M_calc=0;
    phi_G_calc=1;
end

% matrix of the conditional covariance between the past and present:
Cov_X_Y = cov_cond(Cov_X,Cov_XY,Cov_Y); 

% Calculate the normalization factor
% First calculate the differential entropy of each sub-community
num_communities = length(unique(partition));
for i= 1:num_communities
    community_inds = find(partition==i);
    comm_Cov_X{i} = Cov_X(community_inds,community_inds);
    comm_H(i) = entropy_gauss(comm_Cov_X{i});
end
K = (num_communities-1)*min(comm_H); % normalization proposed in Balduzi & Tononi, 2008, PLoS Comp Biol

full_Cov_E = Cov_X_Y; % covariance of the residuals matrix equals the conditional covariance 
%%
if phi_G_calc
    full_A = Cov_XY*inv(Cov_Y); % get the full regression matrix
    output.phi = phi_G_Gauss_AL(Cov_X, full_Cov_E, full_A, partition);
    output.phi_norm = output.phi/K; % normalized phi_G    
elseif phi_M_calc
    % the conditional entropy between the past and present:
    H_cond = entropy_gauss(Cov_X_Y); 
    
    % get the conditional entropies between the past and present of each of
    % your sub-communities
    for i = 1:num_communities
        community_inds = find(partition==i);
        comm_Cov_Y = Cov_Y(community_inds,community_inds);
        comm_Cov_XY = Cov_XY(community_inds,community_inds);
        comm_Cov_X_Y = cov_cond(comm_Cov_X{i},comm_Cov_XY,comm_Cov_Y);
        comm_H_cond(i) = entropy_gauss(comm_Cov_X_Y);
    end
    
    output.phi = sum(comm_H_cond) - H_cond; % phi_M
    output.phi_norm = output.phi/K; % normalized phi_M
end     

end