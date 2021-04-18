function [fingerprint_distance_global, fingerprint_distance_local] = compute_fingerprint_distance(f_phi_1,f_phi_2,n)
% this function calculates the fingerprint distance, area bin-wise (local)
% and spatial graph-spatial graph (global)


    num_bins = numel(f_phi_1);
%     n=2;

    % computing bin-wise f-phi distances
    for i=1:num_bins
      if ~isempty(f_phi_1{1,i}) &&  ~isempty(f_phi_2{1,i})
       bin_distance{i,1} = (abs(f_phi_1{1,i}(:,1) - f_phi_2{1,i}(:,1))).^n;   
       bin_distance{i,2} = f_phi_1{1,i}(:,2);
       bin_distance{i,3} = sum(bin_distance{i,1});
      end 
    end  

    % binwise local distances
    fingerprint_distance_local = cell2mat(bin_distance(:,3));

    % global distances
    fingerprint_distance_global = sum(fingerprint_distance_local.^2);


end

