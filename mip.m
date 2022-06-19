function output = mip(data, phi_measure, tau, extrapolate)

% -----------------------------------------------------------------------
% This code will perform a brute-force search through all possible partitions
% of a system to search for the partition that minimizes normalized
% integrated information. This is EXTRAORDINARILY slow even without
% extrapolating to infinite observations across every partition, so is not
% recommended for any system with more than 8 nodes.
% -----------------------------------------------------------------------

% Generate all possible partitions of the network.
all_partition_inds = SetPartition(num_nodes);
all_partitions = zeros(size(all_partition_inds,1),num_nodes);
for i = 1:size(all_partition_inds,1)
    for j = 1:size(all_partition_inds{i},2)
        all_partitions(i,all_partition_inds{i}{j})=j;
    end
end
all_partitions(1,:)=[];

% Loop through all partitions, computing phi for each,
% and take the smallest result
for i = 1:size(all_partitions,1)
    out = phi_comp(data,all_partitions(i,:),phi_measure,tau,extrapolate);
    min_phi(i) = out.phi;
    min_phi_norm(i) = out.phi_norm;
end
phi_ind = find(phi_all_norm==min(phi_all_norm)); % phi_star across the MIP
output.partition=all_partitions(phi_ind,:);
output.phi=phi_all(phi_ind);
output.phi_norm=phi_all_norm(phi_ind);