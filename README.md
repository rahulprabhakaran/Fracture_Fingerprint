# Fracture_Fingerprint
This repository contains a MATLAB implementation of the fingerprint metric developed by Louf and Barthelemy (2014) for spatial graphs. The fingerprint is a distribution that combines shape factors and block areas that are formed due to the arrangement of spatial graph patterns. Two networks can then be compared by comparing the fingerprint distributions. The Fingerprint_Calculation.m script contains the implementation of the fingerprint distance. Two examples of fracture networks in the form of spatial graphs are included as examples for which the fingerprints are computed and compared. For further information on the fingerprints, please refer to Louf R and Barthelemy M [2014], A typology of street patterns, Journal of The Royal Society Interface, https://doi.org/10.1098/rsif.2014.0924 .

The following data files are available within the repository.
  Fracture Data: g1.mat, xy1.mat, g2.mat, xy2.mat
  Image Data: 20_2_utm299.png, 20_2_utm424.png

We use functions from the following toolboxes which can be downloaded from the MATLAB File Exchange to :
(1) John D'Errico (2021). A suite of minimal bounding objects (https://www.mathworks.com/matlabcentral/fileexchange/34767-a-suite-of-minimal-bounding-objects), MATLAB Central File Exchange. Retrieved July 14, 2020.
(2) Matt J (2021). spatialgraph2D (https://www.mathworks.com/matlabcentral/fileexchange/73630-spatialgraph2d), MATLAB Central File Exchange. Retrieved May 31, 2020.
(3) David Legland (2021). geom2d (https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d), MATLAB Central File Exchange. Retrieved April 04, 2020.
