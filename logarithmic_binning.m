function [freq, bin_values] = logarithmic_binning(data_to_bin, nbins)
% this function bins a given data column vector logarithmically 
% into a specified number of bins, the returned frequency array
% contains the number of observations pertaining to each bin

min_data = min(data_to_bin);
max_data = max(data_to_bin);
data_to_bin_sorted = sort(data_to_bin,'descend');

freq = zeros(nbins,3);
biggest_bin_exp = ceil(log10(max_data));
smallest_bin = 10^biggest_bin_exp/10^(nbins-1);

freq(1,2:3) = [0 smallest_bin];
freq(nbins,2:3) = [10^(biggest_bin_exp-1) 10^biggest_bin_exp];

% assigning number of data points to each bin
for i=1:nbins
   A =  find(data_to_bin>freq(i,2) & data_to_bin<=freq(i,3));
   if isempty(A)
       freq(i,1)=0;
   else
       freq(i,1) = numel(A);
   end
   bin_values{i,1} = data_to_bin(A);
   bin_values{i,2} = A;
   freq(i+1,2) =  freq(i,3);
   freq(i+1,3) =  freq(i,3)*10;
   freq(i,4) = mean([freq(i,2) freq(i,3)]) ;
end
freq(nbins+1,4) = mean([freq(nbins+1,2) freq(nbins+1,3)]) ;
%freq = freq(1:nbins,:);

end

