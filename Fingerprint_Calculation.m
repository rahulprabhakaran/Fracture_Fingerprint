%% Fingerprint Computation
% Rahul Prabhakaran, 2021

% The fingerprint of a spatial graph was introduced by Louf R and
% Barthelemy M [2014], A typology of street patterns, Journal of The Royal 
% Society Interface, https://doi.org/10.1098/rsif.2014.0924 . The same
% paper also presented a fingerprint distance as a means to compare two
% graphs. 

% In this code supplement, we implement the fingerprint distance as
% described by Louf & Barthelemy (2014) for fracture networks that are in
% the form of spatial graphs. 

close all
clear all
format compact
clc

%% loading example fracture data

% load fracture graph and spatial positioning matrix
load('D:\Manuscript_SE\fingerprint_github\g1.mat');
load('D:\Manuscript_SE\fingerprint_github\xy1.mat');

% load image tile extracted from orthomosaic
img = imread('D:\Manuscript_SE\fingerprint_github\20_2_utm424.png');

% computing adjacency matrix for the fracture graph
a = adjacency(g);

% comparing the image tile with the fracture network
figure(1)
subplot(1,2,1)
    imshow(img)
    title('Image tile from Region 1 (1000 x 1000 pixels)')
subplot(1,2,2)
    gplot(a,xy,'k'); pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
    ylim([min(xy(:,2)) max(xy(:,2))]);
    title('Fracture network as spatial graph')
    xlabel('x-coord'); ylabel('y-coord')

%% polygonal face areas computation  

% identifying the polygons from the fracture graph
obj = spatialgraph2D(g,xy(:,1)',xy(:,2)');
[pgon,labelIndices]=polyshape(obj);

blockAreas = zeros(numel(labelIndices),1);
shapeFactor = zeros(numel(labelIndices),1);

% computing block areas and shape factors of each polygon
for i=1:numel(labelIndices)
        disp(i)
        % the third argument is a boolean which uses convex hull to compute the
        % minimum circle of value is 1; if value is zero, then it is not used
        % using '0' can cause errors for some combinations of points, matrices
        % turn out to be not invertible
        [center,radius] = minboundcircle(xy(labelIndices{i},1),xy(labelIndices{i},2),1);
        %figure
        %viscircles(center,radius);
        %hold on
        %scatter(xy(labelIndices{i},1),xy(labelIndices{i},2),'b', 'Filled')
        blockAreas(i,1) = abs(polygonArea(xy(labelIndices{i},1:2)));  
        areaCircumCircle = pi*radius^2;
        shapeFactor(i,1) = blockAreas(i,1)/areaCircumCircle;
        clearvars areaCircumCircle center radius 
end

% converting block areas to cm^2
blockAreas = blockAreas.*1E4;

% binning the block areas into 3 logarithmic bins
[freq, bin_values] = logarithmic_binning(blockAreas, 3);

[bin_1_blocks,~] = find(blockAreas(:,1)>0 & blockAreas(:,1)<= 100);
[bin_2_blocks,~] = find(blockAreas(:,1)>100 & blockAreas(:,1)<= 1000);
[bin_3_blocks,~] = find(blockAreas(:,1)>1000 & blockAreas(:,1)<= 10000);

pgon_bin_1 = pgon(bin_1_blocks);
pgon_bin_2 = pgon(bin_2_blocks);
pgon_bin_3 = pgon(bin_3_blocks);

% plotting fracture network and blocks colored as per block area bins
figure(2)
subplot(1,2,1)
   gplot(a,xy,'k');
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
   title('Fracture network as spatial graph')
   xlabel('x-coord'); ylabel('y-coord')
    
subplot(1,2,2)
   plot(pgon_bin_1,'FaceColor','#0c2c84')  % blue  #0c2c84
   hold on 
   plot(pgon_bin_2,'FaceColor','#fe9929')  % orange #fe9929
   hold on
   plot(pgon_bin_3,'FaceColor','#238b45')  % green #238b45
   axis off
   hold on
   gplot(a,xy,'k')
   title('Block Areas in logarithmic bins')
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
   
   

%% computing fingerprint of a fracture network

fingerprint = [shapeFactor blockAreas];

% logarithmically binning the block areas
% blockAreas = fingerprint(:,2).*10^4;  % converting area from sq.m to sq.cm
% shapeFactor = fingerprint(:,1);

total_area = sum(bin_values{1,1}) + sum(bin_values{2,1}) + sum(bin_values{3,1});

% finding normalized shape factor for block areas in each bin and
% calculating the sum of each normalized curve
sigma_shape_factor =zeros(27,1);

for i=1:numel(bin_values(:,1))
      if ~isempty(bin_values{i,1})  
         bin_values{i,3} = fingerprint(bin_values{i,2},1);   % shapeFactor in col 1
         freq = [];
         freq = [freq; regular_binning(bin_values{i,3},0.04)];
         freq(:,5) = freq(:,1)./numel(bin_values{i,3});
         bin_values{i,4} = [0; freq(:,4); 1];  % interpolated points on x-axis
         bin_values{i,5} = [0; freq(:,5); 0];  % normalized shape factor
         bin_values{i,5} = bin_values{i,5}.*(sum(bin_values{i,1})/total_area*100);
         sigma_shape_factor = sigma_shape_factor + bin_values{i,5};
         clearvars freq
      end
end
sigma_shape_factor = sigma_shape_factor;
[A, ~] = logarithmic_binning(blockAreas, 3);
A = A(1:numel(A(:,1))-1,2:3);

figure(3)
subplot(1,2,1)
   plot(pgon_bin_1,'FaceColor','#0c2c84')  % blue  #0c2c84
   hold on 
   plot(pgon_bin_2,'FaceColor','#fe9929')  % orange #fe9929
   hold on
   plot(pgon_bin_3,'FaceColor','#238b45')  % green #238b45
   axis off
   hold on
   gplot(a,xy,'k')
   title('Block Areas in logarithmic bins')
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
    
subplot(1,2,2)   
    h=area(bin_values{numel(bin_values(:,1)),4},sigma_shape_factor);
    h(1).FaceColor =  [0.85 0.85 0.85]; %[0 0.65 0.65]; 
    h(1).EdgeColor = 'none';
    hold on
    plot(bin_values{3,4},bin_values{3,5},'Color','#00CC00','LineWidth',2);
    hold on
    plot(bin_values{2,4},bin_values{2,5},'Color','#FF9933','LineWidth',2);
    hold on
    plot(bin_values{1,4},bin_values{1,5},'Color','#0080FF','LineWidth',2);
    grid on
    legend('Fingerprint','1000-10000 cm^2','100-1000 cm^2','0-100 cm^2')
    pbaspect([1 1 1])
    xlabel('Phi')
    ylabel('P(Phi)')
    
% computing f_phi pertaining to the fingerprint
    for i=1:numel(bin_values(:,1))
      if ~isempty(bin_values{i,1})  
         bin_values{i,3} = fingerprint(bin_values{i,2},1);   % shapeFactor in col 1
         [bin_values{i,4},bin_values{i,5}] = compute_f_phi(bin_values{i,3},0.01,numel(fingerprint(:,1)));
         f_phi_1{1,i} = [bin_values{i,4} bin_values{i,5}];
      end     
    end      
    
%% loading a different fracture network and computing fingerprint

clearvars a A bin_1_blocks bin_2_blocks bin_3_blocks bin_values blockAreas ...
    fingerprint_distance_global g h i img labelIndices obj pgon pgon_bin_1 ...
    pgon_bin_2 pgon_bin_3 shapeFactor sigma_shape_factor total_area xy

load('D:\Manuscript_SE\fingerprint_github\g2.mat');
load('D:\Manuscript_SE\fingerprint_github\xy2.mat');

img = imread('D:\Manuscript_SE\fingerprint_github\20_2_utm299.png');

% computing adjacency matrix for the fracture graph
a = adjacency(g);

figure(4)
subplot(1,2,1)
    imshow(img)
    title('Image tile from Region 1 (1000 x 1000 pixels)')
subplot(1,2,2)
    gplot(a,xy,'k'); pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
    ylim([min(xy(:,2)) max(xy(:,2))]);
    title('Fracture network as spatial graph')
    xlabel('x-coord'); ylabel('y-coord')

% computing the shape-factor and block areas corresponding to the spatial graph    
obj = spatialgraph2D(g,xy(:,1)',xy(:,2)');
[pgon,labelIndices]=polyshape(obj);

blockAreas = zeros(numel(labelIndices),1);
shapeFactor = zeros(numel(labelIndices),1);

for i=1:numel(labelIndices)
        disp(i)
        % the third argument is a boolean which uses convex hull to compute the
        % minimum circle of value is 1; if value is zero, then it is not used
        % using '0' can cause errors for some combinations of points, matrices
        % turn out to be not invertible
        [center,radius] = minboundcircle(xy(labelIndices{i},1),xy(labelIndices{i},2),1);
        blockAreas(i,1) = abs(polygonArea(xy(labelIndices{i},1:2)));  
        areaCircumCircle = pi*radius^2;
        shapeFactor(i,1) = blockAreas(i,1)/areaCircumCircle;
        clearvars areaCircumCircle center radius 
end

blockAreas = blockAreas.*1E4;

[freq, bin_values] = logarithmic_binning(blockAreas, 3);

[bin_1_blocks,~] = find(blockAreas(:,1)>0 & blockAreas(:,1)<= 100);
[bin_2_blocks,~] = find(blockAreas(:,1)>100 & blockAreas(:,1)<= 1000);
[bin_3_blocks,~] = find(blockAreas(:,1)>1000 & blockAreas(:,1)<= 10000);

pgon_bin_1 = pgon(bin_1_blocks);
pgon_bin_2 = pgon(bin_2_blocks);
pgon_bin_3 = pgon(bin_3_blocks);

figure(5)
subplot(1,2,1)
   gplot(a,xy,'k');
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
   title('Fracture network as spatial graph')
   xlabel('x-coord'); ylabel('y-coord')
    
subplot(1,2,2)
   plot(pgon_bin_1,'FaceColor','#0c2c84')  % blue  #0c2c84
   hold on 
   plot(pgon_bin_2,'FaceColor','#fe9929')  % orange #fe9929
   hold on
   plot(pgon_bin_3,'FaceColor','#238b45')  % green #238b45
   axis off
   hold on
   gplot(a,xy,'k')
   title('Block Areas in logarithmic bins')
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
   
% computing fingerprint   
fingerprint = [shapeFactor blockAreas];

total_area = sum(bin_values{1,1}) + sum(bin_values{2,1}) + sum(bin_values{3,1});

% finding normalized shape factor for block areas in each bin and
% calculating the sum of each normalized curve
sigma_shape_factor =zeros(27,1);

for i=1:numel(bin_values(:,1))
      if ~isempty(bin_values{i,1})  
         bin_values{i,3} = fingerprint(bin_values{i,2},1);   % shapeFactor in col 1
         freq = [];
         freq = [freq; regular_binning(bin_values{i,3},0.04)];
         freq(:,5) = freq(:,1)./numel(bin_values{i,3});

         bin_values{i,4} = [0; freq(:,4); 1];  % interpolated points on x-axis
         bin_values{i,5} = [0; freq(:,5); 0];  % normalized shape factor
         bin_values{i,5} = bin_values{i,5}.*(sum(bin_values{i,1})/total_area*100);
         sigma_shape_factor = sigma_shape_factor + bin_values{i,5};
         clearvars freq
      end
end
sigma_shape_factor = sigma_shape_factor;
[A, ~] = logarithmic_binning(blockAreas, 3);
A = A(1:numel(A(:,1))-1,2:3);

figure(6)
subplot(1,2,1)
   plot(pgon_bin_1,'FaceColor','#0c2c84')  % blue  #0c2c84
   hold on 
   plot(pgon_bin_2,'FaceColor','#fe9929')  % orange #fe9929
   hold on
   plot(pgon_bin_3,'FaceColor','#238b45')  % green #238b45
   axis off
   hold on
   gplot(a,xy,'k')
   title('Block Areas in logarithmic bins')
   pbaspect([1 1 1]); xlim([min(xy(:,1)) max(xy(:,1))]);
   ylim([min(xy(:,2)) max(xy(:,2))]);
    
subplot(1,2,2)   
    h=area(bin_values{numel(bin_values(:,1)),4},sigma_shape_factor);
    h(1).FaceColor =  [0.85 0.85 0.85]; 
    h(1).EdgeColor = 'none';
    hold on
    plot(bin_values{3,4},bin_values{3,5},'Color','#00CC00','LineWidth',2);
    hold on
    plot(bin_values{2,4},bin_values{2,5},'Color','#FF9933','LineWidth',2);
    hold on
    plot(bin_values{1,4},bin_values{1,5},'Color','#0080FF','LineWidth',2);
    grid on
    legend('Fingerprint','1000-10000 cm^2','100-1000 cm^2','0-100 cm^2')
    pbaspect([1 1 1])
    xlabel('Phi')   
    
% computing f_phi pertaining to the fingerprint
    for i=1:numel(bin_values(:,1))
      if ~isempty(bin_values{i,1})  
         bin_values{i,3} = fingerprint(bin_values{i,2},1);   % shapeFactor in col 1
         [bin_values{i,4},bin_values{i,5}] = compute_f_phi(bin_values{i,3},0.01,numel(fingerprint(:,1)));
         f_phi_2{1,i} = [bin_values{i,4} bin_values{i,5}];
      end     
    end  
        

%% comparing two fracture networks using the fingerprint distance

% comparing fingerprints for the same network    
[fingerprint_distance_global, ~] = compute_fingerprint_distance(f_phi_1,f_phi_1,1);  

disp('Fingerprint Distance comparing Network 1 with itself = ')
disp(fingerprint_distance_global)

% comparing fingerprints for the two networks
[fingerprint_distance_global, ~] = compute_fingerprint_distance(f_phi_1,f_phi_2,1);  

disp('Fingerprint Distance comparing Network 1 with Network 2 = ')
disp(fingerprint_distance_global)
