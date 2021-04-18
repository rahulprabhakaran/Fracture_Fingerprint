function [freq] = regular_binning(data_to_bin,increment)
% this function bins normalized data that lies between 0-1
% into regular sized intervals and returns the frequency of
% each bin at a mean point between intervals to enable
% easy plotting


% min_data = min(data_to_bin);
% max_data = max(data_to_bin);
% data_to_bin_sorted = sort(data_to_bin,'descend');


freq(:,2) = 0:increment:1-increment;
freq(:,3) = 0+increment:increment:1;

for i=1:numel(freq(:,1))
   A =  find(data_to_bin>freq(i,2) & data_to_bin<=freq(i,3));
   if isempty(A)
       freq(i,1)=0;
   else
       freq(i,1) = numel(A);
   end   
   freq(i,4)=mean([freq(i,2) freq(i,3)]);
end



end

