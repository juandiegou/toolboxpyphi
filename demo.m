
% -----------------------------------------------------------------------
% This demo will show you how to use this toolbox to calculate integrated
% information from time-series data. 
% -----------------------------------------------------------------------

% LOAD AND CHECK DATA
% load sample data
load('sample_timeseries.mat')

% FIRST MAKE SURE YOUR DATA ARE MULTIVARIATE NORMAL! If they aren't, the
% values you get from anything in this toolbox will *not* be information
% measures. In other words, you can only throw the word "bits" in front of
% the numbers you get if your data are at least approximately
% multivariate Gaussian:

% To check Gaussianity, first plot a multivariate QQ-Plot of your data. If 
% the data are Gaussian, the blue points should mostly fall along the dotted 
% line:
figure;
mv_qqplot(data)

% It's also good to check some histograms of the time-series of
% individual nodes/channels in your data, which should follow a univariate
% normal distribution if your data are multivariate normal:
figure;
hist(data(1,:));

figure; 
hist(data(10,:));

% CALCULATE INTEGRATED INFORMATION: GENERAL HOW-TO
%
% The main function you want to use is phi_comp.m, which calls all the
% other functions you'll need. phi_comp.m takes five inputs, in the
% following order:
%
% 1) Your time-series data (every row is a series of
% observations from one node/channel)
%
% 2) The method you want to use to identify a partition across which to
% measure integrated information. The options are 'MIB' for the minimum
% information bipartition, 'MIP' for the minimum information partition
% (NOT recommended if you want the code to finish running in this
% lifetime), 'SpectClust' for our spectral clustering-based solution,
% 'queyranne' for quickly the Queyranne algorithm (which does a good job
% of minimizing non-normalized integrated information, but here attempts to
% minimize normalized integrated information), 'REMCMC' for a Replica
% Exchange Markov Chain Monte Carlo approach to very, very slowly find
% the MIB (it's slow but tractable approach for large systems), or 'atomic'
% to treat each node as its own community. You can also feed in your own
% community-assignment vector, with just 1s and 2s for a bipartition
%
% 3) The measure of integrated information you want to use. You have two
% options: 1) 'phi_G' for geometric integrated information (which is the
% most theoretically sound option available, and the one we used in the
% paper) or 2) 'phi_M' which is the second most theoretically sound option
% available (its main drawback is that it can exceed the total mutual information 
% in time in a system)
%
% 4) The time-lag across which to calcualte integrated information 
%
% 5) A flag telling the code whether or not to extrapolate integrated
% information to infinite observations across every partition it tries. A 1
% indicates extrapolate, and a 0 indicates don't extrapolate. Extrapolation
% is the more theoretically sound option of the two, but it adds a *lot* to
% your run-time, which is why we didn't use this option in our paper

% CALCULATE INTEGRATED INFORMATION: EXAMPLES

% For the sake of speeding things up, we will not extrapolate to infinite
% observations in the following examples

% Calculate phi_G, non-extrapolated, across the MIB (identified through a
% brute-force search through all possible bipartitions), with a time-lag tau 
% of 3 observations: 
phiG_MIB = phi_comp(data,'MIB','phi_G', 3,0)

% Calculate the same thing, but across the spectral clustering-based
% partition. The result should be identical as it was across the MIB:
phiG_Spect = phi_comp(data,'SpectClust','phi_G', 3,0)

% Now calculate integrated information across the partition identified by
% the Queyranne algorithm. Make sure the 'sfo' folder is in your path.
% This algorithm typically does not find the right answer, since it's
% better suited for minimizing non-normalized integrated information (see
% Kitazono & Oizumi 2017):
phiG_quey = phi_comp(data,'queyranne','phi_G', 3,0)

% Alternately, you can try your own partition:
partition_vector = [1 1 1 2 2 2 2 2 1 1 1 1 1 1];
phiG_part = phi_comp(data, partition_vector, 'phi_G', 3)


% Though we didn't include phi_M in the paper, our code can calculate this
% value as well:
phiM_MIB = phi_comp(data,'MIB','phi_M', 3,0)

% Now go forth and analyze some data!

