% % % %%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
% % % Reconstructs stains from cast-off event to reproduce the motion of cast-off.
% % % 
% % % %D.A. and S.MC acknowledge financial support from the Center for Statistics and Applications in Forensic Evidence (CSAFE), the National Institute of Justice (NIJ), and Iowa State University, as well as pro bono technical and consulting services from Struo LLC, a scientific consulting company based in Ames IA.  D.A. thanks his Iowa State University colleague Xuan Hien Nguyen for pointing the geometry proposition to him.  E.L. would like to thank Independent Forensic Services (IFS) for their assistance and use of their facilities and resources in collecting the non-circular cast-off patterns.
% % % Cast-off Reconstruction has only been tested on Matlab 2019a and GNU Octave 5.2.0. It is recommended that Matlab 2019a, GNU Octave 5.2.0 or newer versions are used to run the included files. Users run the included files at their own risk. Especially, D.A. and S.MC assume no responsibility regarding the results, their interpretation and the consequences thereof.
% % % 
% % % Required Repository Files to run the code:
% % %   - Spatter Measurement Data, e.g. 'Ink_Trial.csv', 'Swineblood_Trial.csv'
% % % - 'Castoff_Reconstruction_DRIVER.m'
% % %  - 'DRIVER.csv' produces 'DRIVER.mat' required for 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_FUNC.m'
% % %  - 'lineSegmentIntersect.m'
% % %  - 'point_to_line.m'
% % %  - 'gauss_distribution.m'
% % %  - 'CircleFitByPratt.m'
% % % - 'Castoff_Reconstruction_POST.m' 
% % %  - 'meshVolume.m' 
% % %  - 'inpolyhedron.m' 
% % %  - 'triangulateFaces.m' 
% % % 
% % % Licenses:
% % % All licenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the repository ‘Licenses’ folder, this was not intentional by the author.
% % % 
% % % Required Installation for GNU Octave Compatibility: (Suggested Installation Order)
% % %  - Updated GNU Octave (written and tested with 5.2.0) (https://www.gnu.org/software/octave/download.html)
% % % 	 - All packages can be updated to the latest version by running:
% % %   		>> pkg update
% % %  - Updated Java (https://www.java.com/en/download/win10.jsp)
% % %  - Psychtoolbox-3 Requirements (install prior to Psychtoolbox download) (http://psychtoolbox.org/download.html#download-problems)
% % % 	 - Psychtoolbox is a free set of Matlab and GNU Octave functions for vision and neuroscience research.
% % % 	 - Subversion 1.7.x command-line client (e.g. Sliksvn) MUST INSTALL COMPLETE SOFTWARE(https://sliksvn.com/download/)
% % % 		 - SlikSVN is a standalone command-line Subversion client for Windows.
% % % 	 - 64-Bit GStreamer-1.16.0 MSVC runtime or later versions MUST INSTALL COMPLETE SOFTWARE(https://gstreamer.freedesktop.org/data/pkg/windows/1.16.0/gstreamer-1.0-msvc-x86_64-1.16.0.msi)		 - GStreamer is a library for constructing graphs of media-handling components.
% % %  - Download Psychtoolbox Toolbox Version 3 (PTB-3) (http://psychtoolbox.org/)
% % % 	 - Unzip download in new folder ‘C:\toolbox’
% % % 	 - Open Octave, copy and paste the following script to the command window to install Psychtoolbox:
% % % >> cd C:\toolbox
% % % >> DownloadPsychtoolbox(‘C:\toolbox’)
% % %  - Follow command window prompts to finish installation.
% % %  - Prior to running ‘Castoff_Recosntruction_DRIVER.m’ or ‘Castoff_Reconstruction_MAIN.m’ the ‘io’ and ‘signal’ packages need to be loaded.
% % % 	 - Load packages on an as needed basis by copying and pasting the following script to the command window:
% % % 		>> pkg load io
% % % 		>> pkg load signal
% % % 	 - Required packages can be found at (https://octave.sourceforge.io/packages.php) for download.
% % % 	 - Automatically load packages at Octave startup with the following:
% % % 		 - Go to (C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup)
% % % 		 - Open ‘octaverc’ with a text editor.
% % % 		 - Copy and paste the following lines to the end of the document:
% % % 		>> pkg load io
% % % 		>> pkg load signal
% % %  - Save file.
% % % 			 - If ‘octaverc’ does not exist, copy, paste, and save the following to a text file and save to the directory (C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup):
% % % 				>> pkg load io
% % % 				>> pkg load signal
% % % 
% % % 
% % % 
% % % Instructions to Run:
% % % 1. Save all included files to the same directory.
% % % 2. Set desired DRIVER (lines 140-259) and MAIN (lines 124-147) User Inputs along with Ad-hoc User defined Variables. (FARO and Hemospat Drivers are provided for trialing code)
% % % 3. Choose desired clustering method ('dwn_samp_stains' and 'opti_space' is the default Clustering Method).
% % % 4. 'Run' (F5) the DRIVER ('Castoff_Reconstruction_DRIVER.m').
% % % 5. 'Run' (F5) the MAIN ('Castoff_Reconstruction_MAIN.m').
% % % 6. MAIN will output the 'Total Elapsed Cluster Analysis Time:', in seconds (s), with the total program run time excluding any variable time from user input. See ‘User Inputs: (DRIVER)’ ‘res’ for estimated runtimes.
% % % 7. MAIN will save the Resultant Variables as a .mat file denoted by ‘Castoff_Reconstruction.mat’ and Output the Cast-off Reconstruction Results as Figures(4+) (figure number is dependent on clustering method). Figure(1) shows the Pratt fir automated reference point. Figure(3) outputs each clustered cast-off reconstructed arc.
% % % 
% % % Figure Displaying:
% % %  - Show and hide legends from resultant figures using:
% % % 		 - legend show
% % % 		 - legend hide
% % % 	 - Rotate, zoom in/out, and pan figure controls are located on the top of the figures to allow for other viewing angles. Also, for all figures besides Figure(1), view(az,el) can be used for precise three-dimensional viewing angles where az is the azimuth angle (in degrees) and el is the polar (or elevation) angle (in degrees).
% % % 
% % % 
% % % User Inputs: (DRIVER)
% % %  - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (lines 147-149)
% % %  - 'lngth','minor' major and minor axis lengths of stains in mm (lines 154-155)
% % %  - 'alpha','gamma' stain impact and directional angles (lines 156-158)
% % %  - 'room', room length (along x-dimension), width (along y-dimension), and height (along z-direction) (x,y,z) in centimeters (line 189)
% % %  - 'xmin','xmax','ymin','ymax','zmin','zmax' minimum and maximum coordinate of possible region of cast-off origin) (lines 183-189)
% % %  - 'res' Spatial resolution of reconstruction (Length of Discretized Uniform Regions of Space Dimensions) (1-15cm is the recommended range) (15cm took ~30 seconds, 10cm took ~1 minute, and 7.5cm took ~1 hour in a large room with several hundred stains and default specifications) (line 193)
% % %  - 'n_b','n_u','n_f','n_d','t_b','t_u','t_f','t_d' surface normal (n) and tangential (t) unit vectors (four surfaces by default are assumed to be perpendicular) (surfaces include: back,upward,front,downward) *** Does not incorporate side surfaces *** (lines 213-220)
% % % - 'inclsurf_[1-4]' Choose '1' to INCLUDE Surface #[1-4] Stains; Choose '0' to EXCLUDE Surface #[1-4]
% % %  - 'InOutTrajectory' *****NOT RECOMMENDED TO USE***** Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
% % %  - 'alpha60more' Choose '1' to Remove impact angles Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees (line 244)
% % %  - 'alphaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Impact Angles (alpha); Choose '0' to NOT Weight Center & Radius Results by Impact Angles (alpha)
% % %  - 'gammaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Glancing Angles (gamma); Choose '0' to NOT Weight Center & Radius Results by Glancing Angles (gamma)
% % %  - 'zetaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Bisector Intersection Angles (zeta); Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
% % %  - 'bsctintweight' *****NOT RECOMMENDED TO USE***** Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta) (line 249)
% % %  - 'InRoom' *****NOT RECOMMENDED TO USE***** Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories
% % % 
% % % User Inputs: (MAIN)
% % %  - 'percentiles' Choose Resultant Product Distribution Percentiles for Plotting Likelihood Regions (this can be change when replotting the final figure (Figure 4) in the Post-processor (line 124)
% % % 
% % % Ad-hoc User Defined Variables:
% % %  - 'stdev' Approximate Standard Deviation of Stain Width Measurement in mm (0.1 is the default value) (DRIVER line 191)
% % %  - 'Spread_Fact_cu' Spreading Factor Uncertainty in Distance between a given Discretized Region in Space and Arc (10cm-30cm is recommended range, 2*res (cm) is the default) (DRIVER line 227)
% % %  - 'Spread_Fact_theta' Spreading Factor Uncertainty in In-plane Angle (Theta) in radians (pi/180 radians (1 degree) is the default) (DRIVER line 228)
% % %  - 'Spread_Fact_upsilon' Spreading Factor Uncertainty in Off-plane Angle (Upsilon) in radians (pi/18-pi/6 radians (10-30 degrees) is recommended pi/9 radians (20 degrees) is the default (DRIVER line 229)
% % %  
% % % Recommended Clustering Methods:
% % %  - 'user_clstr' Set Equal to '1' to Run One Specific Cluster of Three Stains, Set Equal to '0' for another Cluster Option
% % %  - 'stain_cluster' Enter user defined stains indices in an nx3 matrix for n clusters to cluster specific stain combinations
% % %  - 'dwn_samp_stains' Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option (MAIN line 132)
% % %  - 'dwnsamp' Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains (MAIN line 133)
% % %  - 'sampsize' Select Stain Cluster Sample Size by Integer Factor Greater than Three (3) (MAIN line 134)
% % %  - 'overlap' Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired; Set Equal to 'dwnsamp*(sampsize-1)' for maximum overlap (MAIN line 135)
% % %  - 'alpha_clstr' Set Equal to '1' to Cluster by Half Global Alpha Impact Angle, Set Equal to '0' for another Cluster Option
% % %  - 'clstr_alpha' Select Difference in Half Global Alpha Impact Angle for Clustering in radians
% % %  - 'alpha_std' Select Standard Deviation of Difference in Half Global Alpha Impact Angle for Clustering in radians
% % %  - 'dist_clstr' Set Equal to '1' to Cluster by Distance between Stains, Set Equal to '0' for another Cluster Option
% % %  - 'clstr_dist' Select Distance between Stains for Clustering in centimeters
% % %  - 'dist_std' Select Standard Deviation of Distance between Stains for Clustering in centimeters
% % %  - 'opti_space' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1). (MAIN line 146)
% % %  - 'opti_angle' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Angle between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).


%Note to User:
%Script Requires Functions 'Castoff_Reconstruction_DRIVER.m' & 'Castoff_Reconstruction_MAIN.m' & 'Castoff_Reconstruction_FUNC.m' & 'lineSegmentIntersect.m' & 'point_to_line.m' & 'gauss_distribution' & Desired Data, e.g. 'Ink_Trial.csv' or 'Swineblood_Trial.csv' to Run


clear,clc

dbstop if error
warning on verbose

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Enter Actual Center & Radius if Known  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actual_x = 101;%45.72;% %Actual X-Coordinate of Center Location; Enter "NaN" if Unknown
actual_y = 49.21;%48.3;% %Actual X-Coordinate of Center Location; Enter "NaN" if Unknown
actual_z = 95;%89.5;% %Actual Z-Coordinate of Center Location; Enter "NaN" if Unknown
actual_r = 70; %Actual Radius of Cast-off; Enter "NaN" if Unknown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Export Data from Excel Spreadsheet  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drive_data = 'Swineblood_Trial.csv';
[num,txt,raw] = xlsread(drive_data,'A1:O100'); %Read in .CSV data file
drive_datamat = regexprep(drive_data,'.csv','','ignorecase'); %Change .mat name for Saving Results
drive_datamat = strcat(drive_datamat,'_DRIVER.mat'); %Change .mat name for Saving Results

x = num(:,1); %Import X-coordinate Locations in cm
y = num(:,2); %Import Y-coordinate Locations in cm
z = num(:,3); %Import Z-coordinate Locations in cm

x(any(x==0,2))=0.00000000001; %Replace All Zeros with Near Zero Number
y(any(y==0,2))=0.00000000001; %Replace All Zeros with Near Zero Number
z(any(z==0,2))=0.00000000001; %Replace All Zeros with Near Zero Number
lngth = (num(:,5)')'; %Import Major Axis of Stains in mm
minor = (num(:,6)')'; %Import Minor Axis of Stains in mm
alpha = (num(:,7)'*pi/180)'; %Import Alpha Pitching Impact Angle
alpha(any(alpha==0,1))=0.00000000001; %Replace All Zeros with Near Zero Number
gamma = (num(:,9)'*pi/180)'; %Import Gamma Glancing Impact Angle
gamma(any(gamma==0,1))=0.00000000001; %Replace All Zeros with Near Zero Number
gamma(any(gamma>(2*pi),1)) = gamma(any(gamma>(2*pi),1))-2*pi; %Replace Gamma Values Greater than 2*pi Radians with Same Angle within Allotted Zero to 2pi Range
gamma(any(gamma<0,1)) = gamma(any(gamma<0,1))+2*pi; %Replace Gamma Values Less than 0*pi Radians with Same Angle within Allotted Zero to 2pi Range
% for ia = 1:length(gamma)
%     beta(ia) = atan(tan(alpha(ia))/sin(gamma(ia))); %Calculate Beta Yawing Impact Angle
% end
% 
% for i = 1:length(gamma)
%    alpha_p(i) = abs(atan(tan(beta(i))/tan(0.5*pi-gamma(i)))); %Alpha Impact Angle Projected onto XZ-plane for 2D Analysis
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Choose Surfaces to Include  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inclsurf_1 = 1; %Choose '1' to INCLUDE Surface #1 (Back Surface) Stains; Choose '0' to EXCLUDE Surface #1 (Back Surface)
inclsurf_2 = 1; %Choose '1' to INCLUDE Surface #2 (Upward Surface) Stains; Choose '0' to EXCLUDE Surface #2 (Upward Surface)
inclsurf_3 = 1; %Choose '1' to INCLUDE Surface #3 (Front Surface) Stains; Choose '0' to EXCLUDE Surface #3 (Front Surface)
inclsurf_4 = 1; %Choose '1' to INCLUDE Surface #4 (Downward Surface) Stains; Choose '0' to EXCLUDE Surface #4 (Downward Surface)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Define Values  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin = 0.00000000001; %min(x); %Select Minimum X-coordinate for Room Assignment
xmax = 248; %max(x); %Select Maximum X-coordinate for Room Assignment
ymin = 0.00000000001; %min(y); %Select Minimum Y-coordinate for Room Assignment
ymax = 91.4; %max(y); %Select Maximum Y-coordinate for Room Assignment
zmin = 0.00000000001; %min(z); %Select Minimum Z-coordinate for Room Assignment
zmax = 243; %max(z); %Select Maximum Z-coordinate for Room Assignment
aoi = [xmin xmax ymin ymax zmin zmax]; %Room Assignment

stdev = 0.1; %Approximated Standard Deviation of Stain Width Measurement in mm

res = 7.5; %Resolution of Heat Map (Length of Cube Dimensions)

if xmin<xmax
    xcube = xmin:res:xmax; %X-dimension Range of Cube Positions
    ycube = ymin:res:ymax; %Y-dimension Range of Cube Positions
    zcube = zmin:res:zmax; %Z-dimension Range of Cube Positions
else
    xcube = xmin:-res:xmax; %X-dimension Range of Cube Positions
    ycube = ymin:res:ymax; %Y-dimension Range of Cube Positions
    zcube = zmin:res:zmax; %Z-dimension Range of Cube Positions
end

[Xcu,Zcu,Ycu] = meshgrid(xcube,zcube,ycube); %X,Y,Z Meshed Cube Positions (X,Y,Z Minimum of Cube/Bottom,Front,Left Corner)
cu_cx = Xcu+(res/2); %Cube Center X-Coordinate
cu_cy = Ycu+(res/2); %Cube Center Y-Coordinate
cu_cz = Zcu+(res/2); %Cube Center Z-Coordinate
% Xcube = Xcu(:);
% Ycube = Ycu(:);
% Zcube = Zcu(:);

n_b = [-1,0,0]; %Back Surface Normal Vector (Back Surface)
n_u = [0,0,1]; %Upward Surface Normal Vector (Upward Surface)
n_f = [1,0,0]; %Front Surface Normal Vector (Front Surface)
n_d = [0,0,-1]; %Downward Surface Normal Vector (Downward Surface)
t_b = [0,0,-1]; %Back Surface Tangential Vector (Back Surface)
t_u = [1,0,0]; %Upward Surface Tangential Vector (Upward Surface)
t_f = [0,0,-1]; %Front Surface Tangential Vector (Front Surface)
t_d = [1,0,0]; %Downward Surface Tangential Vector (Downward Surface)

room_size = [aoi(1)+0.00000000001 aoi(2) aoi(5)+0.00000000001 aoi(6)]; %Determine Room Size from Room Assignment
max_room_size = max(room_size); %Maximum Room Dimension for Scaling Purposes
min_room_size = min(room_size); %Minimum Room Dimension for Scaling Purposes
Direction = [10*max_room_size 0 10*max_room_size 10*max_room_size]; %Set to [max_room_size 0 max_room_size max_room_size] for Upward cast-off; [max_room_size max_room_size max_room_size 0] for Downward cast-off

Spread_Fact_cu = res*2; %Spreading Factor %Uncertainty in Distance between a given Cube and Arc in centimeters
Spread_Fact_theta = pi/180; %Spreading Factor %Uncertainty in In-plane Angle (Theta) in radians
Spread_Fact_upsilon = 20*pi/180; %Spreading Factor %Uncertainty in Off-plane Angle (Upsilon) in radians

dwnsamp = 4; %Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains
sampsize = 3;%round(0.45*numstains); %Select Stain Cluster Sample Size by Integer Factor Greater than Three (3)
overlap = dwnsamp*(sampsize-1);%sampsize-1; %Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired
A = 2*dwnsamp+(sampsize-2); %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i=1

test = 500; %For Plotting Individual Bisectors with Paired Trajectories *See Lines 466-477

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Settings & Weights  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InOutTrajectory = 0; %Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
% alpha30less = 0; %Choose '1' to Remove Alpha Values Less Than 30 degrees; Choose '0' to Keep All Alpha Values Less Than 30 degrees
alpha60more = 0; %Choose '1' to Remove Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees
gammaweight = 1; %Choose '1' to Weight Center & Radius Results by Glancing Angles (gamma); Choose '0' to NOT Weight Center & Radius Results by Glancing Angles (gamma)
alphaweight = 1; %Choose '1' to Weight Center & Radius Results by Impact Angles (alpha); Choose '0' to NOT Weight Center & Radius Results by Impact Angles (alpha)
% deltaweight = 0; %Choose '1' to Weight Center & Radius Results by Trajectory Intersection Angles (delta); Choose '0' to NOT Weight Center & Radius Results by Trajectory Intersection Angles (delta)
zetaweight = 1; %Choose '1' to Weight Center & Radius Results by Bisector Intersection Angles (zeta); Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
bsctintweight = 0; %Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
InRoom = 0; %Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories

mean_alpha = 50*pi/180; %Mean of Impact Angles for Gaussian Distribution (0-5-10-60-90-120) in RADIANS
stdev_alpha = 5*pi/180; %Standard Deviation of Impact Angles for Gaussian Distribution (180-90) in RADIANS
mean_delta = 90*pi/180; %Mean of Trajectory Intersection Angles for Gaussian Distribution (0-5-10-60-90-120) in RADIANS
stdev_delta = 15*pi/180; %Standard Deviation of Trajectory Intersection Angles for Gaussian Distribution (180-90) in RADIANS
mean_zeta = 90*pi/180; %Mean of Bisector Trajectory Intersection Angles for Gaussian Distribution (0-5-10-60-90-120) in RADIANS
stdev_zeta = 30*pi/180; %Standard Deviation of Bisector Trajectory Intersection Angles for Gaussian Distribution (180-90) in RADIANS
mean_bsctint = 0; %Mean of Distance between Trajectory Intersection & Bisector Intersection for Gaussian Distribution (0-5-10-60-90-120) in RADIANS
stdev_bsctint = max_room_size; %Standard Deviation of Distance between Trajectory Intersection & Bisector Intersection for Gaussian Distribution (180-90) in RADIANS

if alpha60more == 1
    x(any(alpha>(60*pi/180),1)) = [];
    y(any(alpha>(60*pi/180),1)) = [];
    z(any(alpha>(60*pi/180),1)) = [];
    gamma(any(alpha>(60*pi/180),1)) = [];
    beta(any(alpha>(60*pi/180),1)) = [];
    alpha_p(any(alpha>(60*pi/180),1)) = [];
    alpha(any(alpha>(60*pi/180),1)) = [];
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Define Room Lengths in Centimeters & Resolution Cubes  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = abs(xmax-xmin); %Room Length in X-direction (cm)
Ly = abs(ymax-ymin); %Room Length in Y-direction (cm)
Lz = abs(zmax-zmin); %Room Length in Z-direction (cm)
Nx = ceil(Lx/res); %Room Length in X-direction (Resolution Cubes)
Ny = ceil(Ly/res); %Room Length in Y-direction (Resolution Cubes)
Nz = ceil(Lz/res); %Room Length in Z-direction (Resolution Cubes)
isocubes = zeros(Nz,Nx,Ny); %Preallocate Cubes
isocubes = isocubes(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Reorganize Intersections from Stains to Surfaces  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Surface Intersection Coordinates
surface_B1 = [room_size(2) room_size(3) room_size(2) room_size(4)]; %Surface #1 (Front) in [X1 Y1 X2 Y2] coordinates
surface_U2 = [room_size(2) room_size(4) room_size(1) room_size(4)]; %Surface #2 (Upward) in [X1 Y1 X2 Y2] coordinates
surface_F3 = [room_size(1) room_size(4) room_size(1) room_size(3)]; %Surface #3 (Back) in [X1 Y1 X2 Y2] coordinates
surface_D4 = [room_size(1) room_size(3) room_size(2) room_size(3)]; %Surface #4 (Downward) in [X1 Y1 X2 Y2] coordinates

surf1 = (xmin-2<=x & x<=xmin+2)'; %Surface Criteria
surf2 = (zmax-2<=z & z<=zmax+2)'; %Surface Criteria
surf3 = (xmax-2<=x & x<=xmax+2)'; %Surface Criteria
surf4 = (zmin-2<=z & z<=zmin+2)'; %Surface Criteria

for ia = 1:length(gamma)
    beta(ia,:) = atan(tan(alpha(ia))/sin(gamma(ia))); %Calculate Beta Yawing Impact Angle
end

for i = 1:length(gamma)
   alpha_p(i,:) = abs(atan(tan(beta(i))/tan(gamma(i)))); %Alpha Impact Angle Projected onto XZ-plane for 2D Analysis
end

if any(surf1(:))==0
    inclsurf_1 = 0;
end
if any(surf2(:))==0
    inclsurf_2 = 0;
end
if any(surf3(:))==0
    inclsurf_3 = 0;
end
if any(surf4(:))==0
    inclsurf_4 = 0;
end

%Remove Excluded Surface Stains as User Previously Defined
if inclsurf_1 ~= 1
    surf1 = zeros(size(surf1)); %Remove Surface #1 Stains
end
if inclsurf_2 ~= 1
    surf2 = zeros(size(surf2)); %Remove Surface #2 Stains
end
if inclsurf_3 ~= 1
    surf3 = zeros(size(surf3)); %Remove Surface #3 Stains
end
if inclsurf_4 ~= 1
    surf4 = zeros(size(surf4)); %Remove Surface #4 Stains
end

% surfaceF1 = [nonzeros(surf1.*x'), nonzeros(surf1.*z'), nonzeros(surf1.*alpha_p), nonzeros(surf1.*gamma)]; %Apply Surface Criteria to Surfaces
% surfaceU2 = [nonzeros(surf2.*x'), nonzeros(surf2.*z'), nonzeros(surf2.*alpha_p), nonzeros(surf2.*gamma)]; %Apply Surface Criteria to Surfaces
% surfaceB3 = [nonzeros(surf3.*x'), nonzeros(surf3.*z'), nonzeros(surf3.*alpha_p), nonzeros(surf3.*gamma)]; %Apply Surface Criteria to Surfaces
% surfaceD4 = [nonzeros(surf4.*x'), nonzeros(surf4.*z'), nonzeros(surf4.*alpha_p), nonzeros(surf4.*gamma)]; %Apply Surface Criteria to Surfaces
% surfaceB1 = [surf1'.*x, surf1'.*y, surf1'.*z, surf1'.*alpha_p', surf1'.*gamma', surf1'.*beta', surf1'.*alpha', surf1'.*minor', surf1'.*lngth']; %Apply Surface Criteria to Surfaces
% surfaceU2 = [surf2'.*x, surf2'.*y, surf2'.*z, surf2'.*alpha_p', surf2'.*gamma', surf2'.*beta', surf2'.*alpha', surf2'.*minor', surf2'.*lngth']; %Apply Surface Criteria to Surfaces
% surfaceF3 = [surf3'.*x, surf3'.*y, surf3'.*z, surf3'.*alpha_p', surf3'.*gamma', surf3'.*beta', surf3'.*alpha', surf3'.*minor', surf3'.*lngth']; %Apply Surface Criteria to Surfaces
% surfaceD4 = [surf4'.*x, surf4'.*y, surf4'.*z, surf4'.*alpha_p', surf4'.*gamma', surf4'.*beta', surf4'.*alpha', surf4'.*minor', surf4'.*lngth']; %Apply Surface Criteria to Surfaces
% surfaceB1 = sortrows(surfaceB1,3); %Reorder Surface Vectors into Chronological Order
% surfaceU2 = sortrows(surfaceU2,1); %surfaceU2 = flip(surfaceU2,1); %Reorder Surface Vectors into Chronological Order
% surfaceF3 = sortrows(surfaceF3,3); surfaceF3 = flip(surfaceF3,1); %Reorder Surface Vectors into Chronological Order
% surfaceD4 = sortrows(surfaceD4,1); surfaceD4 = flip(surfaceD4,1); %Reorder Surface Vectors into Chronological Order
% surfaceB1(isnan(surfaceB1)) = 0; %Set NaNs Equal to Zero
% surfaceU2(isnan(surfaceU2)) = 0; %Set NaNs Equal to Zero
% surfaceF3(isnan(surfaceF3)) = 0; %Set NaNs Equal to Zero
% surfaceD4(isnan(surfaceD4)) = 0; %Set NaNs Equal to Zero
% Xs1 = surfaceB1(:,1); %Select Surface Endpoint
% Xs2 = surfaceU2(:,1); %Select Surface Endpoint
% Xs3 = surfaceF3(:,1); %Select Surface Endpoint
% Xs4 = surfaceD4(:,1); %Select Surface Endpoint
% Ys1 = surfaceB1(:,2); %Select Surface Endpoint
% Ys2 = surfaceU2(:,2); %Select Surface Endpoint
% Ys3 = surfaceF3(:,2); %Select Surface Endpoint
% Ys4 = surfaceD4(:,2); %Select Surface Endpoint
% Zs1 = surfaceB1(:,3); %Select Surface Endpoint
% Zs2 = surfaceU2(:,3); %Select Surface Endpoint
% Zs3 = surfaceF3(:,3); %Select Surface Endpoint
% Zs4 = surfaceD4(:,3); %Select Surface Endpoint
% alpha_pxz1 = surfaceB1(:,4); %Select Surface Projected Alpha Angle
% alpha_pxz2 = surfaceU2(:,4); %Select Surface Projected Alpha Angle
% alpha_pxz3 = surfaceF3(:,4); %Select Surface Projected Alpha Angle
% alpha_pxz4 = surfaceD4(:,4); %Select Surface Projected Alpha Angle
% gamma_1 = surfaceB1(:,5); %Select Surface Gamma Angle
% gamma_2 = surfaceU2(:,5); %Select Surface Gamma Angle
% gamma_3 = surfaceF3(:,5); %Select Surface Gamma Angle
% gamma_4 = surfaceD4(:,5); %Select Surface Gamma Angle
% beta_1 = surfaceB1(:,6); %Select Surface Beta Angle
% beta_2 = surfaceU2(:,6); %Select Surface Beta Angle
% beta_3 = surfaceF3(:,6); %Select Surface Beta Angle
% beta_4 = surfaceD4(:,6); %Select Surface Beta Angle
% alpha_1 = surfaceB1(:,7); %Select Surface Alpha Angle
% alpha_2 = surfaceU2(:,7); %Select Surface Alpha Angle
% alpha_3 = surfaceF3(:,7); %Select Surface Alpha Angle
% alpha_4 = surfaceD4(:,7); %Select Surface Alpha Angle
% minor_1 = surfaceB1(:,8); %Select Surface Minor Axis
% minor_2 = surfaceU2(:,8); %Select Surface Minor Axis
% minor_3 = surfaceF3(:,8); %Select Surface Minor Axis
% minor_4 = surfaceD4(:,8); %Select Surface Minor Axis
% lngth_1 = surfaceB1(:,9); %Select Surface Major Axis
% lngth_2 = surfaceU2(:,9); %Select Surface Major Axis
% lngth_3 = surfaceF3(:,9); %Select Surface Major Axis
% lngth_4 = surfaceD4(:,9); %Select Surface Major Axis
% 
% Xs = [Xs1, Xs2, Xs3, Xs4]; %Compile Endpoints
% Ys = [Ys1, Ys2, Ys3, Ys4]; %Compile Endpoints
% Zs = [Zs1, Zs2, Zs3, Zs4]; %Compile Endpoints
% alpha_p = [alpha_pxz1, alpha_pxz2, alpha_pxz3, alpha_pxz4]; %Compile Angles
% gamma = [gamma_1, gamma_2, gamma_3, gamma_4]; %Compile Angles
% beta = [beta_1, beta_2, beta_3, beta_4]; %Compile Angles
% alpha = [alpha_1, alpha_2, alpha_3, alpha_4]; %Compile Angles
% minor = [minor_1, minor_2, minor_3, minor_4]; %Compile Minor Axes
% lngth = [lngth_1, lngth_2, lngth_3, lngth_4]; %Compile Minor Axes

Xs = x;
Ys = y;
Zs = z;

int_s1 = sum(surf1);
int_s2 = sum(surf2);
int_s3 = sum(surf3);
int_s4 = sum(surf4);

% int_s1 = size(nonzeros(surfaceB1(:,1)),1); %Number of Trajectory Intersections with Surface #1
% int_s2 = size(nonzeros(surfaceU2(:,1)),1); %Number of Trajectory Intersections with Surface #2
% int_s3 = size(nonzeros(surfaceF3(:,1)),1); %Number of Trajectory Intersections with Surface #3
% int_s4 = size(nonzeros(surfaceD4(:,1)),1); %Number of Trajectory Intersections with Surface #4

inclsurf = sum([inclsurf_1,inclsurf_2,inclsurf_3,inclsurf_4].*[1,2,3,4],1); inclsurf = sum(inclsurf,2);
surfmat = [[1,2,3,4];[2,3,4,1];[1,2,3,4];[1,2,3,4];[4,1,2,3];[1,2,3,4];[3,4,1,2];[3,4,1,2];[2,3,4,1];[1,2,3,4]];

% Xs = nonzeros([Xs(:,surfmat(inclsurf,1));Xs(:,surfmat(inclsurf,2));Xs(:,surfmat(inclsurf,3));Xs(:,surfmat(inclsurf,4))]);
% Ys = nonzeros([Ys(:,surfmat(inclsurf,1));Ys(:,surfmat(inclsurf,2));Ys(:,surfmat(inclsurf,3));Ys(:,surfmat(inclsurf,4))]);
% Zs = nonzeros([Zs(:,surfmat(inclsurf,1));Zs(:,surfmat(inclsurf,2));Zs(:,surfmat(inclsurf,3));Zs(:,surfmat(inclsurf,4))]);
% alpha_p = nonzeros([alpha_p(:,surfmat(inclsurf,1));alpha_p(:,surfmat(inclsurf,2));alpha_p(:,surfmat(inclsurf,3));alpha_p(:,surfmat(inclsurf,4))]);
% gamma = nonzeros([gamma(:,surfmat(inclsurf,1));gamma(:,surfmat(inclsurf,2));gamma(:,surfmat(inclsurf,3));gamma(:,surfmat(inclsurf,4))]);
% beta = nonzeros([beta(:,surfmat(inclsurf,1));beta(:,surfmat(inclsurf,2));beta(:,surfmat(inclsurf,3));beta(:,surfmat(inclsurf,4))]);
% alpha = nonzeros([alpha(:,surfmat(inclsurf,1));alpha(:,surfmat(inclsurf,2));alpha(:,surfmat(inclsurf,3));alpha(:,surfmat(inclsurf,4))]);
% minor = nonzeros([minor(:,surfmat(inclsurf,1));minor(:,surfmat(inclsurf,2));minor(:,surfmat(inclsurf,3));minor(:,surfmat(inclsurf,4))]);
% lngth = nonzeros([lngth(:,surfmat(inclsurf,1));lngth(:,surfmat(inclsurf,2));lngth(:,surfmat(inclsurf,3));lngth(:,surfmat(inclsurf,4))]);

numstains = size(Xs,1); %Revaluating Variables

%Surface Vector [1n; 2n; 3n; 4n;]
for j = 1:numstains
    if room_size(1)-2<=Xs(j) && Xs(j)<=room_size(1)+2
        face(j) = 1; %Stains Contained on Surface #1
    elseif room_size(2)-2<=Xs(j) && Xs(j)<=room_size(2)+2
        face(j) = 3; %Stains Contained on Surface #3
    elseif room_size(3)-2<=Zs(j) && Zs(j)<=room_size(3)+2
        face(j)=4; %Stains Contained on Surface #4
    elseif room_size(4)-2<=Zs(j) && Zs(j)<=room_size(4)+2
        face(j)=2; %Stains Contained on Surface #2
    end
end

face = face'; %Row to Column Vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  alpha_p to alpha_pg  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Using gamma to Identify Upward/Downward Stains, Denoted as +/-
for j = 1:numstains
        if 0<=gamma(j) && gamma(j)<0.5*pi
            alpha_p(j) = -1*(alpha_p(j)); %Define Projected Alpha Direction by Gamma and Surface          
        elseif 0.5*pi<=gamma(j) && gamma(j)<0.5*3*pi
            alpha_p(j) = alpha_p(j);     %Define Projected Alpha Direction by Gamma and Surface
        elseif 0.5*3*pi<=gamma(j) && gamma(j)<=2*pi
            alpha_p(j) = -1*(alpha_p(j)); %Define Projected Alpha Direction by Gamma and Surface
        else
            alpha_p(j) = NaN; %Define Projected Alpha Direction by Gamma and Surface
        end   
end

rmvnan = double(~isnan(alpha_p)); %Projected Alpha Corrections
rmvnan(rmvnan==0) = NaN; %Prjected Alpha Corrections

%Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
for k = 1:numstains
    if alpha_p<0
        if 0<=k & k<=int_s1
            alpha_pg(k,:) = pi+abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif int_s1<k & k<=(int_s1+int_s2)
            alpha_pg(k,:) = 0.5*pi-abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif (int_s1+int_s2)<k & k<=(int_s1+int_s2+int_s3)
            alpha_pg(k,:) = abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif (int_s1+int_s2+int_s3)<k & k<=(int_s1+int_s2+int_s3+int_s4)
            alpha_pg(k,:) = 0.5*3*pi-abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        end
    else
        if 0<=k & k<=int_s1
            alpha_pg(k,:) = 2*pi-abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif int_s1<k & k<=(int_s1+int_s2)
            alpha_pg(k,:) = 0.5*3*pi+abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif (int_s1+int_s2)<k & k<=(int_s1+int_s2+int_s3)
            alpha_pg(k,:) = pi-abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        elseif (int_s1+int_s2+int_s3)<k & k<=(int_s1+int_s2+int_s3+int_s4)
            alpha_pg(k,:) = 0.5*pi+abs(alpha_p(k)); %Calculate Projected Alpha with Respect to Gravity from Projected Alpha and Surface
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Matlab Structuring  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matlab Structuring
for jj=1:numstains
    castoff(jj).name=strcat('Stain_#',int2str(jj)); %Unique Identifier for Each Stain Vector
    castoff(jj).S=face(jj); %Surface vector
    castoff(jj).X_cm = Xs(jj); %X-Coordinate on Surface Vector
    castoff(jj).Y_cm = Ys(jj); %Y-Coordinate on Surface Vector
    castoff(jj).Z_cm = Zs(jj); %Z-Coordinate on Surface Vector
    castoff(jj).ALPHA_PG_rad = alpha_pg(jj); %Projected Alpha_gravity (Angle between Point on Cast-off Circle, Intersection of Trajectory & Surface, and Plumb Reference Vector
    castoff(jj).ALPHA_rad = alpha(jj); %Alpha Pitching Impact Angle
    castoff(jj).WIDTH_mm = minor(jj); %Minor Axis of Stains in mm
    castoff(jj).LENGTH_mm = lngth(jj); %Major Axis of Stains in mm
    castoff(jj).CENTER_cm = [actual_x,actual_y,actual_z]; %Armature Center Coordinates
    castoff(jj).RADIUS_cm = actual_r; %Armature Radius
    castoff(jj).DIST_cm = sqrt((Xs(jj)-actual_x)^2+(Ys(jj)-actual_y)^2+(Zs(jj)-actual_z)^2); %Distance between Stain and Armature Center
end

numstains = length(x); %Number of Stains in .CSV  

save(drive_datamat);

fprintf('Script ran at %s\n', datestr(now,'HH:MM:SS.FFF')); %End Run Designator with Time
toc
