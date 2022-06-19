function output = spectral_partition(data, phi_measure, tau, extrapolate)

% -----------------------------------------------------------------------
% This function implements the spectral clustering-based solution to
% finding the MIB described in the paper. It takes the correlation matrix
% of your time-series data, and transforms it across a sweep of two
% parameters: the steepness or a power law adjacency transformation and
% binarization threshold. Each such transformation yields one guess for
% the "functional network" that describes your data (a concept borrowed
% from neural connectomics). Spectral clustering is then applied to each
% transformed correlation matrix, and integrated information is calculated
% across each resulting bipartition. The function's guess for the MIB of
% your data is the bipartition (among these) that minimizes normalized
% integrated information.
% -----------------------------------------------------------------------



% Parameters to sweep through:

% A vector of exponents for the power adjacency transformation of the
% correlation matrix. The exponenets are logarithmically spaced from 0
% to 10
pow_vec = [0 logspace((log(1)/log(10)),(log(20)/log(10)),10)];

% A vector of percentiles at which to binarize each transformed correlation 
% matrix. For faster (but sometimes less accurate) run-time, this vector
% can be changed to include fewer values (e.g. quant_vec=0:.1:.9)
quant_vec=0:.005:.99;

% Get the correlation matrix of your time-series data
corr_vals = corr(data');


if extrapolate 
    % (The code looks a bit different in the extrapolated and non-extrapolated
    % cases, because you'll need to calculate different covariance matrices 
    % for each sub-sample of your data)
    
    % Create cells that will store all candidate partitions, values of
    % integrated information, and values of normalized integrated
    % information:
    all_spect_partitions = cell(1,length(pow_vec));
    all_phi = cell(1,length(pow_vec));
    all_phi_norm = cell(1,length(pow_vec));
    
    % Sweep through power law adjacency transformations:
    for pow = 1:length(pow_vec)
        if pow_vec(pow)==0
            adj=abs(corr_vals); 
        else
            adj = ((corr_vals+1)./2).^pow_vec(pow); % transformation
        end
        % Sweep through cut-offs for binarization:
        for i = 1:length(quant_vec)
            spect_A = adj;
            if quant_vec(i)>0 % no binarization if 0
                cutoff = quantile(adj(:),quant_vec(i));
                spect_A(adj<cutoff)=0;         
            end
            % Apply spectral clustering to your transformed correlation
            % matrix, splitting your system into 2 communities, and repeat
            % k-means clustering 3 times:
            VV= GCSpectralClust2(spect_A,2,3); 
            spect_partitions(i,:)=VV(:,2)';
            
            % Calculate integrated information across the spectral
            % partition:
            out = extrapolated_phi_calc(data,spect_partitions(i,:),phi_measure,tau);
            min_phi(i) = out.phi; % integrated information
            min_phi_norm(i) = out.phi_norm; % normalized integrated information 
        end
        all_spect_partitions{pow} = spect_partitions; % store all candidate partitions
        all_phi{pow} = min_phi; % store all values of integrated information
        all_phi_norm{pow} = min_phi_norm; % store all values of normalized integrated information
    end
else 
    % (In the non-extrapolated case, we only have to set up our covariance
    % matrices once, and feed them into phi_calc.m as opposed to
    % phi_comp.m, as we do above for the extrapolated case. This saves a
    % lot of computation time)
    
    nobs = size(data,2); % number of observations/time steps in your data
    % Define the time-ranges for the past and present of your system:
    timerange1 = 1:1:nobs-tau; % past
    timerange2 = 1+tau:1:nobs; % present
    % Covariance matrix for the present of your system:
    Cov_X = data(:,timerange2)*data(:,timerange2)'/(nobs-tau-1);
    % Covariance matrix for the past of your system:
    Cov_Y = data(:,timerange1)*data(:,timerange1)'/(nobs-tau-1); 
    % Cross-covariance between the past and present of your system:
    Cov_XY = (data(:,timerange2)*data(:,timerange1)')/(nobs-tau-1); 
    
    
    % Create cells that will store all candidate partitions, values of
    % integrated information, and values of normalized integrated
    % information:
    all_spect_partitions = cell(1,length(pow_vec));
    all_phi = cell(1,length(pow_vec));
    all_phi_norm = cell(1,length(pow_vec));
    
    % Sweep through power law adjacency transformations:
    for pow = 1:length(pow_vec)
        % Apply the power adjacency transformation to the correlation
        % matrix
        if pow_vec(pow)==0
            adj=abs(corr_vals); 
        else
            adj = ((corr_vals+1)./2).^pow_vec(pow); 
        end
        
        % Sweep through cut-offs for binarization:
        for i = 1:length(quant_vec)
            spect_A = adj;
            if quant_vec(i)>0 % no binarization if 0
                cutoff = quantile(adj(:),quant_vec(i));
                spect_A(adj<cutoff)=0;             
            end
            
            % Apply spectral clustering to your transformed correlation
            % matrix, splitting your system into 2 communities, and repeat
            % k-means clustering 3 times:
            VV= GCSpectralClust2(spect_A,2,3);
            spect_partitions(i,:)=VV(:,2)';
            
            % Calculate integrated information across the spectral
            % partition:
            out = phi_calc(phi_measure,Cov_X,Cov_Y,Cov_XY,spect_partitions(i,:));
            min_phi(i) = out.phi; % integrated information
            min_phi_norm(i) = out.phi_norm; % normalized integrated information 
        end
        all_spect_partitions{pow} = spect_partitions; % store all partitions
        all_phi{pow} = min_phi; % store integrated information values
        all_phi_norm{pow} = min_phi_norm; % store normalized integrated information values
    end
end

all_spect_partitions = cell2mat(all_spect_partitions');
all_phi = cell2mat(all_phi);
all_phi_norm = cell2mat(all_phi_norm);

% Find the partition that minimized normalized integrated information:
min_phi_ind = find(all_phi_norm==min(all_phi_norm));

% Spectral clustering will often produce the exact MIB for a range of
% parameters in this sweep, so just pick the index for the first time the
% algorithm landed on the partition that minimized integrated information: 
if length(min_phi_ind)>1
    min_phi_ind=min_phi_ind(1);
end

% Output:
output.partition = all_spect_partitions(min_phi_ind,:); % best guess for the MIB
output.phi = all_phi(min_phi_ind); % integrated information across the best guess for the MIB
output.phi_norm = all_phi_norm(min_phi_ind); % normalized integrated information across the best guess for the MIB