function [A,B] = compute_f_phi(shapefactor_data,increment,num_blocks)
% computes ratio of number of blocks with a shape factor within a required 
% interval specified by the increment argument to total number of blocks; 
% the returned value is a two column matrix in which first column has the
% ratio and second column is interpolated centre of intervals for plotting


    N1N2(:,1) = 0:increment:1-increment;
    N1N2(:,2) = 0+increment:increment:1;


    for i=1:numel(N1N2(:,1)) 
     f_phi(i,1)= numel(find(shapefactor_data>N1N2(i,1) & shapefactor_data<=N1N2(i,2)))/num_blocks;
     f_phi(i,2)= mean([N1N2(i,1) N1N2(i,2)]); 
    end 

    f_phi_plot = [0 0];
    f_phi_plot = [f_phi_plot; f_phi];
    
    f_phi_plot = [f_phi_plot; [0 1]];
    A = f_phi_plot(:,1);
    B = f_phi_plot(:,2);
end

