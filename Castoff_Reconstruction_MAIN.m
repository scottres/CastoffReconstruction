% % % %%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
% % % Reconstructs stains from cast-off event to reproduce the motion of cast-off.
% % %
% % % Last Updated 07/27/2020
% % % 
% % % %D.A. and S.MC acknowledge financial support from the Center for Statistics and Applications in Forensic Evidence (CSAFE), the National Institute of Justice (NIJ), and Iowa State University, as well as pro bono technical and consulting services from Struo LLC, a scientific consulting company based in Ames IA.  D.A. thanks his Iowa State University colleague Xuan Hien Nguyen for pointing the geometry proposition to him.  E.L. would like to thank Independent Forensic Services (IFS) for their assistance and use of their facilities and resources in collecting the non-circular cast-off patterns.
% % % Cast-off Reconstruction has only been tested on Matlab 2019a and GNU Octave 5.2.0. It is recommended that Matlab 2019a, GNU Octave 5.2.0 or newer versions are used to run the included files. Users run the included files at their own risk. Especially, D.A. and S.MC assume no responsibility regarding the results, their interpretation and the consequences thereof.
% % % 
% % % Required Repository Files to run the code:
% % %  - Spatter Measurement Data, e.g. 'Ink_Trial.csv', 'Swineblood_Trial.csv'
% % %  - 'Castoff_Reconstruction_DRIVER.m'
% % %  - 'DRIVER.csv' produces 'DRIVER.mat' required for 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_FUNC.m'
% % %  - 'lineSegmentIntersect.m'
% % %  - 'point_to_line.m'
% % %  - 'gauss_distribution.m'
% % %  - 'CircleFitByPratt.m'
% % %  - 'Castoff_Reconstruction_POST.m' 
% % %  - 'meshVolume.m' 
% % %  - 'inpolyhedron.m' 
% % %  - 'triangulateFaces.m' 
% % % 
% % % Licenses:
% % % All licenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the repository ˜Licenses folder, this was not intentional by the author.
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
% % % 	 - Unzip download in new folder ˜C:\toolbox
% % % 	 - Open Octave, copy and paste the following script to the command window to install Psychtoolbox:
% % % >> cd C:\toolbox
% % % >> DownloadPsychtoolbox(˜C:\toolbox)
% % %  - Follow command window prompts to finish installation.
% % %  - Prior to running ˜Castoff_Recosntruction_DRIVER.m or ˜Castoff_Reconstruction_MAIN.m the ˜io and ˜signal packages need to be loaded.
% % % 	 - Load packages on an as needed basis by copying and pasting the following script to the command window:
% % % 		>> pkg load io
% % % 		>> pkg load signal
% % % 	 - Required packages can be found at (https://octave.sourceforge.io/packages.php) for download.
% % % 	 - Automatically load packages at Octave startup with the following:
% % % 		 - Go to (C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup)
% % % 		 - Open ˜octaverc with a text editor.
% % % 		 - Copy and paste the following lines to the end of the document:
% % % 		>> pkg load io
% % % 		>> pkg load signal
% % %  - Save file.
% % % 			 - If ˜octaverc does not exist, copy, paste, and save the following to a text file and save to the directory (C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup):
% % % 				>> pkg load io
% % % 				>> pkg load signal
% % % 
% % % Instructions to Run:
% % % 1. Save all included files to the same directory.
% % % 2. Set desired DRIVER (lines 152-270) and MAIN (lines 134-158) User Inputs along with Ad-hoc User defined Variables. (FARO and Hemospat Drivers are provided for trialing code)
% % % 3. Choose desired clustering method ('dwn_samp_stains' and 'opti_space' is the default Clustering Method).
% % % 4. 'Run' (F5) the DRIVER ('Castoff_Reconstruction_DRIVER.m').
% % % 5. 'Run' (F5) the MAIN ('Castoff_Reconstruction_MAIN.m').
% % % 6. MAIN will output the 'Total Elapsed Cluster Analysis Time:', in seconds (s), with the total program run time excluding any variable time from user input. See ˜User Inputs: (DRIVER) ˜res for estimated runtimes.
% % % 7. MAIN will save the Resultant Variables as a .mat file denoted by ‘Castoff_Reconstruction.mat’ and Output the Cast-off Reconstruction Results with warnings as a txt-file denoted by 'Castoff_Reconstruction_OUTPUT.txt' displayed in Figures(4+) (figure number is dependent on clustering method). Figure(1) shows the Pratt fit automated reference point. Figure(3) outputs each clustered cast-off reconstructed arc.
% % % 
% % % Figure Displaying:
% % %  - Show and hide legends from resultant figures using:
% % % 		 - legend show
% % % 		 - legend hide
% % % 	 - Rotate, zoom in/out, and pan figure controls are located on the top of the figures to allow for other viewing angles. Also, for all figures besides Figure(1), view(az,el) can be used for precise three-dimensional viewing angles where az is the azimuth angle (in degrees) and el is the polar (or elevation) angle (in degrees).
% % % 
% % % User Inputs: (DRIVER)
% % %  - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (lines 158-160)
% % %  - 'lngth','minor' major and minor axis lengths of stains in mm (lines 165-166)
% % %  - 'alpha','gamma' stain impact and directional angles (lines 167-169)
% % %  - 'room', room length (along x-dimension), width (along y-dimension), and height (along z-direction) (x,y,z) in centimeters (line 200)
% % %  - 'xmin','xmax','ymin','ymax','zmin','zmax' minimum and maximum coordinate of possible region of cast-off origin) (lines 194-200)
% % %  - 'res' Spatial resolution of reconstruction (Length of Discretized Uniform Regions of Space Dimensions) (1-15cm is the recommended range) (15cm took ~30 seconds, 10cm took ~1 minute, and 7.5cm took ~1 hour in a large room with several hundred stains and default specifications) (line 204)
% % %  - 'n_b','n_u','n_f','n_d','t_b','t_u','t_f','t_d' surface normal (n) and tangential (t) unit vectors (four surfaces by default are assumed to be perpendicular) (surfaces include: back,upward,front,downward) *** Does not incorporate side surfaces *** (lines 224-231)
% % % - 'inclsurf_[1-4]' Choose '1' to INCLUDE Surface #[1-4] Stains; Choose '0' to EXCLUDE Surface #[1-4]
% % %  - 'InOutTrajectory' *****NOT RECOMMENDED TO USE***** Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
% % %  - 'alpha60more' Choose '1' to Remove impact angles Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees (line 255)
% % %  - 'alphaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Impact Angles (alpha); Choose '0' to NOT Weight Center & Radius Results by Impact Angles (alpha)
% % %  - 'gammaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Glancing Angles (gamma); Choose '0' to NOT Weight Center & Radius Results by Glancing Angles (gamma)
% % %  - 'zetaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Bisector Intersection Angles (zeta); Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
% % %  - 'bsctintweight' *****NOT RECOMMENDED TO USE***** Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta) (line 260)
% % %  - 'InRoom' *****NOT RECOMMENDED TO USE***** Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories
% % % 
% % % User Inputs: (MAIN)
% % %  - 'percentiles' Choose Resultant Product Distribution Percentiles for Plotting Likelihood Regions (this can be change when replotting the final figure (Figure 4) in the Post-processor (line 133)
% % % 
% % % Ad-hoc User Defined Variables:
% % %  - 'stdev' Approximate Standard Deviation of Stain Width Measurement in mm (0.1 is the default value) (DRIVER line 202)
% % %  - 'Spread_Fact_cu' Spreading Factor Uncertainty in Distance between a given Discretized Region in Space and Arc (10cm-30cm is recommended range, 2*res (cm) is the default) (DRIVER line 238)
% % %  - 'Spread_Fact_theta' Spreading Factor Uncertainty in In-plane Angle (Theta) in radians (pi/180 radians (1 degree) is the default) (DRIVER line 239)
% % %  - 'Spread_Fact_upsilon' Spreading Factor Uncertainty in Off-plane Angle (Upsilon) in radians (pi/18-pi/6 radians (10-30 degrees) is recommended pi/9 radians (20 degrees) is the default (DRIVER line 240)
% % %  
% % % Recommended Clustering Methods:
% % %  - 'user_clstr' Set Equal to '1' to Run One Specific Cluster of Three Stains, Set Equal to '0' for another Cluster Option
% % %  - 'stain_cluster' Enter user defined stains indices in an nx3 matrix for n clusters to cluster specific stain combinations
% % %  - 'dwn_samp_stains' Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option (MAIN line 143)
% % %  - 'dwnsamp' Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains (MAIN line 144)
% % %  - 'sampsize' Select Stain Cluster Sample Size by Integer Factor Greater than Three (3) (MAIN line 145)
% % %  - 'overlap' Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired; Set Equal to 'dwnsamp*(sampsize-1)' for maximum overlap (MAIN line 146)
% % %  - 'alpha_clstr' Set Equal to '1' to Cluster by Half Global Alpha Impact Angle, Set Equal to '0' for another Cluster Option
% % %  - 'clstr_alpha' Select Difference in Half Global Alpha Impact Angle for Clustering in radians
% % %  - 'alpha_std' Select Standard Deviation of Difference in Half Global Alpha Impact Angle for Clustering in radians
% % %  - 'dist_clstr' Set Equal to '1' to Cluster by Distance between Stains, Set Equal to '0' for another Cluster Option
% % %  - 'clstr_dist' Select Distance between Stains for Clustering in centimeters
% % %  - 'dist_std' Select Standard Deviation of Distance between Stains for Clustering in centimeters
% % %  - 'opti_space' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1). (MAIN line 157)
% % %  - 'opti_angle' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Angle between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).

clear,clc,close all

test_interpreter=["Is the interpreter Matlab or " "Octave?"];
if (size(test_interpreter) == [1 36])
  fprintf('Octave interpreter detected.') % Octave is probably running.
  pkg load io
  pkg load signal 
elseif (size(test_interpreter) == [1 2])
  fprintf('Matlab interpreter detected.') %MATLAB is probably running
else
  fprintf('It looks like the interpreter is neither Matlab nor Octave')
end

dbstop if error
warning on verbose

tic()

percentiles = [0.95 0.80 0.65]; %Choose Resultant Product Distribution Percentiles for Plotting (this can be change when replotting the final figure (Figure 4)
data = 'stain Scott3_DRIVER.mat';
load(data); %'Castoff_Reconstruction_DRIVER_SWINEBLOOD.mat' OR 'Castoff_Reconstruction_DRIVER_LISCIO_TRIAL_2_UPDATED.mat'
datamat = regexprep(data,'_DRIVER','','ignorecase'); %Change .mat name for Saving Results
dlmwrite(regexprep(datamat,'.mat','_OUTPUT.txt'),[]); %Create Output .txt-file

user_clstr = 0; %Set Equal to '1' to Run One Specific Cluster of Three Stains, Set Equal to '0' for another Cluster Option
stain_cluster = [52 30 11]; %Enter user defined stains indices in an nx3 matrix for n clusters to cluster specific stain combinations

dwn_samp_stains = 1; %Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option
dwnsamp = 1;%[1:15]'; %Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains
sampsize = 3;%round(0.45*numstains); %Select Stain Cluster Sample Size by Integer Factor Greater than Three (3)
overlap = dwnsamp*(sampsize-1);%sampsize-1; %Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired
A = 2*dwnsamp+(sampsize-2); %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i=1

alpha_clstr = 0; %Set Equal to '1' to Cluster by Half Global Alpha Impact Angle, Set Equal to '0' for another Cluster Option
clstr_alpha = [1*pi/180:2*pi/180:115*pi/180]';%[20*pi/180; 30*pi/180]; %Select Difference in Half Global Alpha Impact Angle for Clustering in radians
alpha_std = 1*pi/180; %Select Standard Deviation of Difference in Half Global Alpha Impact Angle for Clustering in radians

dist_clstr = 0; %Set Equal to '1' to Cluster by Distance between Stains, Set Equal to '0' for another Cluster Option
clstr_dist = 25;%[2:4:100]'; %Select Distance between Stains for Clustering in centimeters
dist_std = 5; %Select Standard Deviation of Distance between Stains for Clustering in centimeters

opti_space = 1; %Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).
opti_angle = 0; %Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Angle between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).

SF_cu_range = [0,inf]; %Distance Spreading Factor Allowable Range
SF_theta_range = [0,2*pi]; %In-plane Angle Spreading Factor Allowable Range
SF_upsilon_range = [0,2*pi]; %Off-plane Angle Spreading Factor Allowable Range
dalpha_range = [0,0.5*pi]; %Alpha Impact Angle Range
dgamma_range = [0,2*pi]; %Gamma Directional Angle Range
res_range = [0,15]; %Spatial Region Resolution Allowable Range
user_default = 0; %Ideal Clustering Value
dwn_samp_default = 1; %Ideal Clustering Value
dwnsamp_default = 1; %Ideal Clustering Value
sampsize_default = 3; %Ideal Clustering Value
overlap_default = dwnsamp*(sampsize-1); %Ideal Clustering Value
alpha_default = 0; %Ideal Clustering Value
dist_default = 0; %Ideal Clustering Value
opti_default = 1; %Ideal Clustering Value

%Input Variables and Warnings for OUTPUT
fileID = fopen(regexprep(datamat,'.mat','_OUTPUT.txt'),'w'); %Open OUTPUT
fprintf(fileID, '%s\r\n\n', regexprep(datamat,'.mat','_OUTPUT.txt'),''); %Output Title to OUTPUT
fprintf(fileID, '%s\r\n\n', 'INPUTS:'); %Output Inputs Header
if Spread_Fact_cu >= SF_cu_range(1) && Spread_Fact_cu < SF_cu_range(2)
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_cu =','',num2str(Spread_Fact_cu))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "Spread_Fact_cu" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_cu =','',num2str(Spread_Fact_cu),' (ideal "Spread_Fact_cu" Range =','',num2str(SF_cu_range))); %Inputs are within Predefined Range
end
if Spread_Fact_theta >= SF_theta_range(1) && Spread_Fact_theta < SF_theta_range(2)
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_theta =','',num2str(Spread_Fact_theta))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "Spread_Fact_theta" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_theta =','',num2str(Spread_Fact_theta),' (ideal "Spread_Fact_theta" Range =','',num2str(SF_theta_range))); %Inputs are within Predefined Range
end
if Spread_Fact_upsilon >= SF_upsilon_range(1) && Spread_Fact_upsilon < SF_upsilon_range(2)
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_upsilon =','',num2str(Spread_Fact_upsilon))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "Spread_Fact_upsilon" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('Spread_Fact_upsilon =','',num2str(Spread_Fact_upsilon),' (ideal "Spread_Fact_upsilon" Range =','',num2str(SF_upsilon_range))); %Inputs are within Predefined Range
end
if res >= res_range(1) && res < res_range(2)
  fprintf(fileID, '%s\r\n', strcat('res =','',num2str(res))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "res" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('res =','',num2str(res),' (ideal "res" Range =','',num2str(res_range))); %Inputs are within Predefined Range
end
if user_clstr == user_default
  fprintf(fileID, '%s\r\n', strcat('user_clstr =','',num2str(user_clstr))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "user_clstr" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('user_clstr =','',num2str(user_clstr),' (ideal "user_clstr" value =','',num2str(user_default),')')); %Inputs are within Predefined Range
end
if dwn_samp_stains == dwn_samp_default
  fprintf(fileID, '%s\r\n', strcat('dwn_samp_stains =','',num2str(dwn_samp_stains))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "dwn_samp_stains" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('dwn_samp_stains =','',num2str(dwn_samp_stains),' (ideal "dwn_samp_stains" value =','',num2str(dwn_samp_default),')')); %Inputs are within Predefined Range
end
if dwnsamp == dwnsamp_default
  fprintf(fileID, '%s\r\n', strcat('dwnsamp =','',num2str(dwnsamp))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "dwnsamp" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('dwnsamp =','',num2str(dwnsamp),' (ideal "dwnsamp" value =','',num2str(dwnsamp_default),')')); %Inputs are within Predefined Range
end
if sampsize == sampsize_default
  fprintf(fileID, '%s\r\n', strcat('sampsize =','',num2str(sampsize))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "sampsize" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('sampsize =','',num2str(sampsize),' (ideal "sampsize" value =','',num2str(sampsize_default),')')); %Inputs are within Predefined Range
end
if overlap == overlap_default
  fprintf(fileID, '%s\r\n', strcat('overlap =','',num2str(overlap))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "overlap" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('overlap =','',num2str(overlap),' (ideal "overlap" value =','',num2str(overlap_default),')')); %Inputs are within Predefined Range
end
if alpha_clstr == alpha_default
  fprintf(fileID, '%s\r\n', strcat('alpha_clstr =','',num2str(alpha_clstr))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "alpha_clstr" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('alpha_clstr =','',num2str(alpha_clstr),' (ideal "alpha_clstr" value =','',num2str(alpha_default),')')); %Inputs are within Predefined Range
end
if dist_clstr == dist_default
  fprintf(fileID, '%s\r\n', strcat('dist_clstr =','',num2str(dist_clstr))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "dist_clstr" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('dist_clstr =','',num2str(dist_clstr),' (ideal "dist_clstr" value =','',num2str(dist_default),')')); %Inputs are within Predefined Range
end
if opti_space == opti_default
  fprintf(fileID, '%s\r\n', strcat('opti_space =','',num2str(opti_space))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "opti_space" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('opti_space =','',num2str(opti_space),' (ideal "opti_space" value =','',num2str(opti_default),')')); %Inputs are within Predefined Range
end
fclose(fileID); %Close OUTPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Define Surface Normal Vectors  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Surface Normal Vector Directions (Normal Directed Outside Room)
nxt_b = cross(n_b,t_b); %Normal x Tangential Vector (Back Surface)
nxt_u = cross(n_u,t_u); %Normal x Tangential Vector (Upward Surface)
nxt_f = cross(n_f,t_f); %Normal x Tangential Vector (Front Surface)
nxt_d = cross(n_d,t_d); %Normal x Tangential Vector (Downward Surface)

en_b = [n_b(1)*surf1' n_b(2)*surf1' n_b(3)*surf1']; %Normal Vectors on Back Surface
en_u = [n_u(1)*surf2' n_u(2)*surf2' n_u(3)*surf2']; %Normal Vectors on Upward Surface
en_f = [n_f(1)*surf3' n_f(2)*surf3' n_f(3)*surf3']; %Normal Vectors on Front Surface
en_d = [n_d(1)*surf4' n_d(2)*surf4' n_d(3)*surf4']; %Normal Vectors on Downward Surface

et_b = [t_b(1)*surf1' t_b(2)*surf1' t_b(3)*surf1']; %Tangential Vectors on Back Surface
et_u = [t_u(1)*surf2' t_u(2)*surf2' t_u(3)*surf2']; %Tangential Vectors on Upward Surface
et_f = [t_f(1)*surf3' t_f(2)*surf3' t_f(3)*surf3']; %Tangential Vectors on Front Surface
et_d = [t_d(1)*surf4' t_d(2)*surf4' t_d(3)*surf4']; %Tangential Vectors on Downward Surface

enxt_b = [nxt_b(1)*surf1' nxt_b(2)*surf1' nxt_b(3)*surf1']; %Normal x Tangential Vectors on Back Surface
enxt_u = [nxt_u(1)*surf2' nxt_u(2)*surf2' nxt_u(3)*surf2']; %Normal x Tangential Vectors on Upward Surface
enxt_f = [nxt_f(1)*surf3' nxt_f(2)*surf3' nxt_f(3)*surf3']; %Normal x Tangential Vectors on Front Surface
enxt_d = [nxt_d(1)*surf4' nxt_d(2)*surf4' nxt_d(3)*surf4']; %Normal x Tangential Vectors on Downward Surface

en = cat(3,en_b,en_u,en_f,en_d); %Concatenate Normal Vectors
en1 = [en(:,:,surfmat(inclsurf,1))]; %Reorder Surface Normal Vectors for Clustering
en2 = [en(:,:,surfmat(inclsurf,2))]; %Reorder Surface Normal Vectors for Clustering
en3 = [en(:,:,surfmat(inclsurf,3))]; %Reorder Surface Normal Vectors for Clustering
en4 = [en(:,:,surfmat(inclsurf,4))]; %Reorder Surface Normal Vectors for Clustering
en1(all(en1==0,2),:)=[]; %Remove Rows of All Zero Elements
en2(all(en2==0,2),:)=[]; %Remove Rows of All Zero Elements
en3(all(en3==0,2),:)=[]; %Remove Rows of All Zero Elements
en4(all(en4==0,2),:)=[]; %Remove Rows of All Zero Elements

et = cat(3,et_b,et_u,et_f,et_d); %Concatenate Tangential Vectors
et1 = [et(:,:,surfmat(inclsurf,1))]; %Reorder Surface Tangential Vectors for Clustering
et2 = [et(:,:,surfmat(inclsurf,2))]; %Reorder Surface Tangential Vectors for Clustering
et3 = [et(:,:,surfmat(inclsurf,3))]; %Reorder Surface Tangential Vectors for Clustering
et4 = [et(:,:,surfmat(inclsurf,4))]; %Reorder Surface Tangential Vectors for Clustering
et1(all(et1==0,2),:)=[]; %Remove Rows of All Zero Elements
et2(all(et2==0,2),:)=[]; %Remove Rows of All Zero Elements
et3(all(et3==0,2),:)=[]; %Remove Rows of All Zero Elements
et4(all(et4==0,2),:)=[]; %Remove Rows of All Zero Elements

enxt = cat(3,enxt_b,enxt_u,enxt_f,enxt_d); %Concatenate Normal x Tangential Vectors
enxt1 = [enxt(:,:,surfmat(inclsurf,1))]; %Reorder Surface Normal x Tangential Vectors for Clustering
enxt2 = [enxt(:,:,surfmat(inclsurf,2))]; %Reorder Surface Normal x Tangential Vectors for Clustering
enxt3 = [enxt(:,:,surfmat(inclsurf,3))]; %Reorder Surface Normal x Tangential Vectors for Clustering
enxt4 = [enxt(:,:,surfmat(inclsurf,4))]; %Reorder Surface Normal x Tangential Vectors for Clustering
enxt1(all(enxt1==0,2),:)=[]; %Remove Rows of All Zero Elements
enxt2(all(enxt2==0,2),:)=[]; %Remove Rows of All Zero Elements
enxt3(all(enxt3==0,2),:)=[]; %Remove Rows of All Zero Elements
enxt4(all(enxt4==0,2),:)=[]; %Remove Rows of All Zero Elements

e_n = sum(en,3); %Normal Vectors
e_t = sum(et,3); %Tangential Vectors
e_nxt = sum(enxt,3); %Normal x Tangential Vectors

%Surface Gamma Vector Directions (For Vertical Surfaces: Gamma = 0 In Direction of Gravity, For Horizontal Surfaces: Gamma = 0 In +X Direction)
%Find Vector Components for Each Stain
v_n = sin(alpha); %Normal Velocity Vector
v_nxt = (v_n./tan(beta)); %Normal x Tangential Velocity Vector
v_t = (v_nxt./tan(gamma)); %Tangential Velocity Vector
v_n = [e_n(:,1).*v_n e_n(:,2).*v_n e_n(:,3).*v_n]; %Velocity Vector in Direction of Normal Vector
v_nxt = [e_nxt(:,1).*v_nxt e_nxt(:,2).*v_nxt e_nxt(:,3).*v_nxt]; %Velocity Vector in Direction of Normal x Tangential Vector
v_t = [e_t(:,1).*v_t e_t(:,2).*v_t e_t(:,3).*v_t]; %Velocity Vector in Direction of Tangential Vector
v = [v_n+v_t+v_nxt]; %Complete Velocity Vector

x_orig = x; %Save Original Stain X-coordinates
y_orig = y; %Save Original Stain Y-coordinates
z_orig = z; %Save Original Stain Z-coordinates
v_orig = v; %Save Original Stain V Vectors
numstains_orig = numstains; %Save Original Number of Stains

%Calculate Average Distance between Known Cast-off Center Location and Adjacent Stains
xi1 = 1; %Starting Stain Index
xij(1,:) = xi1; %Stain Iterative Index
opti_dist = mean(mean(sqrt(sum(diff([Xs Ys Zs],1,1).^2,2)),1),2); %Average Distance between Adjacent Stains (if stains are not inputted in consecutive order this will not be accurate and less stains will be selected for analysis)
% opti_dist = 1; %opti_dist*0.25; %Manually Selected Distance between Adjacent Stains (if stains are not inputted in consecutive order this will not be accurate and less stains will be selected for analysis)

%Select Stains Equally Spaced by 'opti_dist' Distance between Stains
while xi1 < numstains; % && (xi1+25) <= (numstains-1);
    if (xi1+floor(numstains/3))-1 < numstains; %Iterating through the next floor(numstains/3)-1 consecutive stains to prevent Skipping over sections of stains (when cast-off spatter patterns pass by one another or intersect)
        xnums = xi1+floor(numstains/3)-1;
    else
        xnums = numstains-1;
    end
    for jdi = xi1:xnums; %(numstains-1); %(xi1+25);
        dist(jdi,:) = sqrt(sum([(Xs(xi1)-Xs(jdi+1)),(Ys(xi1)-Ys(jdi+1)),(Zs(xi1)-Zs(jdi+1))].^2,2)); %Distance between the next 25 consecutive stains
        close_dist(jdi,:) = abs(opti_dist-dist(jdi,:)); %Choosing stain closest to 'opti_dist' distance from indexed stain
    end;
    close_dist(close_dist == 0) = NaN; %Avoid Selecting Errored Stains
    [Mi,Ii] = min(close_dist); %Determine Stain Index of Selected Stain
    xij(xi1+1,:) = Ii+1; %Saving Stain Index
    xi1 = Ii+1; %Selecting Next Stain Index for Next Iteration
    clearvars dist close_dist %Remove Variables to Avoid Overwriting Errors
end

xij = nonzeros(xij); %Remove Skipped Stains

if length(xij) == 2;
  xij = [xij(1); round((xij(2)-xij(1))/2); xij(2)];
end

%Calculate Average Angle between Known Cast-off Center Location and
%Adjacent Stains *** ONLY APPLICABLE WHEN KNOWN CAST-OFF MOTION IS KNOWN ***
if opti_angle == 1;
    yi1 = 1; %Starting Stain Index
    yij(1,:) = yi1; %Stain Iterative Index

    for ci = 1:(numstains-1);
        abi(ci,:) = [(actual_x-Xs(ci)),(actual_y-Ys(ci)),(actual_z-Zs(ci))]; %Vector between Known Cast-off Motion Center and Indexed Stain
        cbi(ci,:) = [(actual_x-Xs(ci+1)),(actual_y-Ys(ci+1)),(actual_z-Zs(ci+1))]; %Vector between Known Cast-off Motion Center and Adjacent Stain Index
        abci(ci,:) = acos(dot(abi(ci,:),cbi(ci,:))/(sqrt(sum((abi(ci,:)).^2,2))*sqrt(sum((cbi(ci,:)).^2,2)))); %Angle between Vectors 'abi' and 'cbi'
    end

    opti_ang = mean(abci); %Average Angle between Adjacent Stains (if stains are not inputted in consecutive order this will not be accurate and less stains will be selected for analysis)
    
%Select Stains Equally Spaced by 'opti_ang' Angle between Stains
    while yi1 < numstains
        for cj = yi1:(numstains-1);
            abj(cj,:) = [(actual_x-Xs(yi1)),(actual_y-Ys(yi1)),(actual_z-Zs(yi1))]; %Vector between Known Cast-off Motion Center and Indexed Stain for Each Remaining Stain
            cbj(cj,:) = [(actual_x-Xs(cj+1)),(actual_y-Ys(cj+1)),(actual_z-Zs(cj+1))]; %Vector between Known Cast-off Motion Center and Adjacent Stain Index for Each Remaining Stain
            abcj(cj,:) = acos(dot(abj(cj,:),cbj(cj,:))/(sqrt(sum((abj(cj,:)).^2,2))*sqrt(sum((cbj(cj,:)).^2,2)))); %Angle between Vectors 'abi' and 'cbi' for Each Remaining Stain
            close_ang(cj,:) = abs(opti_ang-abcj(cj,:)); %Difference between Angles and 'opti_ang' Angle
        end
        close_ang(close_ang == 0) = NaN; %Avoid Selecting Errored Stains
        [Mj,Ij] = min(close_ang); %Determine Stain Index of Selected Stain
        yij(yi1+1,:) = Ij+1; %Saving Stain Index
        yi1 = Ij+1; %Stain Iterative Index
        clearvars abj cbj abcj close_ang %Remove Variables to Avoid Overwriting Errors
    end

    yij = nonzeros(yij); %Remove Skipped Stains
end

if opti_space == 1 %Only Select Stains Equally Spaced 'opti_dist' Distance Apart
    Xs = Xs(xij,:);
    Ys = Ys(xij,:);
    Zs = Zs(xij,:);
    gamma = gamma(xij,:);
    beta = beta(xij,:);
    alpha_p = alpha_p(xij,:);
    alpha = alpha(xij,:);
    face = face(xij,:);
    alpha_pg = alpha_pg(xij,:);
    minor = minor(xij,:);
    lngth = lngth(xij,:);
    surf1 = surf1(xij);
    surf2 = surf2(xij);
    surf3 = surf3(xij);
    surf4 = surf4(xij);
    en_b = en_b(xij,:);
    en_u = en_u(xij,:);
    en_f = en_f(xij,:);
    en_d = en_d(xij,:);
    et_b = et_b(xij,:);
    et_u = et_u(xij,:);
    et_f = et_f(xij,:);
    et_d = et_d(xij,:);
    enxt_b = enxt_b(xij,:);
    enxt_u = enxt_u(xij,:);
    enxt_f = enxt_f(xij,:);
    enxt_d = enxt_d(xij,:);
    e_n = e_n(xij,:);
    e_t = e_t(xij,:);
    e_nxt = e_nxt(xij,:);
    v_n = v_n(xij,:);
    v_t = v_t(xij,:);
    v_nxt = v_nxt(xij,:);
    v = v(xij,:);
end

if opti_angle == 1 %Only Select Stains Equally Spaced 'opti_ang' Angle Apart
    Xs = Xs(yij,:);
    Ys = Ys(yij,:);
    Zs = Zs(yij,:);
    gamma = gamma(yij,:);
    beta = beta(yij,:);
    alpha_p = alpha_p(yij,:);
    alpha = alpha(yij,:);
    face = face(yij,:);
    alpha_pg = alpha_pg(yij,:);
    minor = minor(yij,:);
    lngth = lngth(yij,:);
    surf1 = surf1(yij);
    surf2 = surf2(yij);
    surf3 = surf3(yij);
    surf4 = surf4(yij);
    en_b = en_b(yij,:);
    en_u = en_u(yij,:);
    en_f = en_f(yij,:);
    en_d = en_d(yij,:);
    et_b = et_b(yij,:);
    et_u = et_u(yij,:);
    et_f = et_f(yij,:);
    et_d = et_d(yij,:);
    enxt_b = enxt_b(yij,:);
    enxt_u = enxt_u(yij,:);
    enxt_f = enxt_f(yij,:);
    enxt_d = enxt_d(yij,:);
    e_n = e_n(yij,:);
    e_t = e_t(yij,:);
    e_nxt = e_nxt(yij,:);
    v_n = v_n(yij,:);
    v_t = v_t(yij,:);
    v_nxt = v_nxt(yij,:);
    v = v(yij,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Automatically Determine Bisector Reference Point  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine Reference Point within Region of Cast-off Origin
xzi = [Xs Ys Zs];%[x y z]+max_room_size*v; %Initial Trajectory Endpoint
xzf = [Xs Ys Zs]-max_room_size*v; %Final Trajectory Endpoint
xrefi = combvec(xzi(:,1)',xzi(:,1)')'; %Compiling All Possible Trajectory Final X-Coordinate Combinations and Purmutations
zrefi = combvec(xzi(:,3)',xzi(:,3)')'; %Compiling All Possible Trajectory Final Z-Coordinate Combinations and Purmutations
xreff = combvec(xzf(:,1)',xzf(:,1)')'; %Compiling All Possible Trajectory Final X-Coordinate Combinations and Purmutations
zreff = combvec(xzf(:,3)',xzf(:,3)')'; %Compiling All Possible Trajectory Final Z-Coordinate Combinations and Purmutations
% zrefi(any(diff(sort(xrefi,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
% xreff(any(diff(sort(xrefi,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
% zreff(any(diff(sort(xrefi,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
% xrefi(any(diff(sort(xrefi,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations

xzrefi = [xrefi(:,1) zrefi(:,1) xreff(:,1) zreff(:,1)]; %Compile Trajectories
xzreff = [xrefi(:,2) zrefi(:,2) xreff(:,2) zreff(:,2)]; %Compile Adjacent Trajectories
xzrefi(isinf(xzrefi) | isnan(xzrefi)) = NaN; %Remove Infinities
xzreff(isinf(xzreff) | isnan(xzreff)) = NaN; %Remove Infinities
xzreff(any(isnan(xzrefi)==1,2),:) = []; %Remove NaNs
xzrefi(any(isnan(xzrefi)==1,2),:) = []; %Remove NaNs
xzrefi(any(isnan(xzreff)==1,2),:) = []; %Remove NaNs
xzreff(any(isnan(xzreff)==1,2),:) = []; %Remove NaNs

%Iterating through all Trajectories to Determine Trajectory Intersections
for xzk =1:length(xzrefi)
    lineSegmentIntersect(xzrefi(xzk,:),xzreff(xzk,:)); %Function Determining Line Segement Intersections: [Xc1 Zc1 Xs1 Zs1] intersecting [Xc2 Zc2 Xs2 Zs2]
    xref(xzk,:) = ans.intMatrixX; %X-Coordinate Trajectory Intersection between Adjacent & Trajectory traj_shift Intersection
    zref(xzk,:) = ans.intMatrixY; %Z-Coordinate Trajectory Intersection between Adjacent & Trajectory traj_shift Intersection
end

xref = nonzeros(xref); %Remove Zeros
zref = nonzeros(zref); %Remove Zeros
xref = xref(~isnan(xref)); %Remove NaNs
zref = zref(~isnan(zref)); %Remove NaNs

xzref = [xref zref]; %Save Circle Center X and Z-coordinates
XZref = CircleFitByPratt(xzref); %Calculate Geometric Fitted Circlular Arc from Stain Trajectory Intersections using the Pratt Method
x_ref_xz = XZref(1); %Save Arc Center X-coordinate
y_ref_yz = mean(Ys); %Select Arc Center Y-coordinate as Average Stain Y-coordinate
z_ref_xz = XZref(2); %Save Arc Center Z-coordinate

Ref = [x_ref_xz y_ref_yz z_ref_xz]; %Compiled User Defined Refernce Point
x_ref = Ref(1); %X-coordinate of User Defined Refernce Point
y_ref = Ref(2); %Y-coordinate of User Defined Refernce Point
z_ref = Ref(3); %Z-coordinate of User Defined Refernce Point

%Plot Room Dimensions, Stains, Stain Trajectories, and Automatically Generated Reference Point
figure(1)
hold on 
grid on
axis equal
h1_front_xz = plot([aoi(1) aoi(1)], [aoi(5) aoi(6)],'Color','c','LineWidth',5); %Plot Front Surface Dimensions
h1_downward_xz = plot([aoi(1) aoi(2)], [aoi(5) aoi(5)],'Color','g','LineWidth',4); %Plot Downward Surface Dimensions
h1_back_xz = plot([aoi(2) aoi(2)], [aoi(5) aoi(6)],'Color','b','LineWidth',3); %Plot Back Surface Dimensions
h1_upward_xz = plot([aoi(1) aoi(2)], [aoi(6) aoi(6)],'Color','y','LineWidth',2); %Plot Downward Surface Dimensions
h2_xz = plot(Xs,Zs, '.','MarkerSize',max_room_size*0.1,'Color','r','LineWidth',2); %Plot XYZ Stains
title({'Cast-off Center & Radius of Origin',''})
xlabel(['X-Axis (cm)'])
ylabel(['Z-Axis (cm)'])

sample1 = numstains; %Save 'numstains' for Debugging
numstains = length(Xs); %Recalculate 'numstains'

for ref1 = 1:numstains
    h3_xz(ref1) = plot([(Xs(ref1)-1000*v(ref1,1)) (Xs(ref1))],[(Zs(ref1)-1000*v(ref1,3)) (Zs(ref1))]); %Plot Stain Trajectories
end
axis square;
xlim([min_room_size-100,max_room_size+100]);
ylim([min_room_size-100,max_room_size+100]);
zlim([min_room_size-100,max_room_size+100]);
% plot(Xref,Zref,'c.','MarkerSize',20); %Plot Automatically Generated Reference Point from Geometric Circle Fit Method
h_guess = plot(XZref(1),XZref(2),'g.','MarkerSize',20);
h_xz = [h1_front_xz h1_downward_xz h1_back_xz h1_upward_xz h2_xz h3_xz(1) h_guess]; % h9(sample2) h10 %Plot Automatically Generated Reference Point by Pratt Method
legend(h_xz, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Trajectories','Automated Reference Point by Pratt Method', 'Location', 'northeastoutside')
chi = linspace(0, 2*pi, 25); %Angle Vector for Plotting Known Cast-off Motion

% ************************************************
% % ***** Use when Cast-off Motion is known *****
% x_fin = (center(1) + radius.*cos(chi))'; %Resultant X-Coordinate of Cast-off Circle
% y_fin = (center(2)*ones(size(x_fin)));
% z_fin = (center(3) + radius.*sin(chi))'; %Resultant Z-Coordinate of Cast-off Circle

% x_actual = actual_x + actual_r.*cos(chi); %Resultant X-Coordinate of Cast-off Circle
% y_actual = actual_y*ones(size(x_actual)); %Resultant Y-Coordinate of Cast-off Circle
% z_actual = actual_z + actual_r.*sin(chi); %Resultant Z-Coordinate of Cast-off Circle

% h15 = plot(x_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
% h16 = plot(actual_x,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
% ************************************************

figure(2); 
hold on; 
h1_front = plot3([aoi(1) aoi(1) aoi(1) aoi(1) aoi(1)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','c','LineWidth',5); 
h1_downward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(5) aoi(5) aoi(5) aoi(5) aoi(5)],'Color','g','LineWidth',4); 
h1_back = plot3([aoi(2) aoi(2) aoi(2) aoi(2) aoi(2)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','b','LineWidth',3); 
h1_upward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)],[aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(6) aoi(6) aoi(6) aoi(6) aoi(6)],'Color','y','LineWidth',2); 
h2 = plot3(x_orig,y_orig,z_orig, '.','MarkerSize',25,'Color','r'); 
xlabel(['X-Axis (cm)']); 
ylabel(['Y-Axis (cm)']); 
zlabel(['Z-Axis (cm)']); 
for ref1 = 1:numstains_orig;    
  h11(ref1) = plot3([(x_orig(ref1)-10000*v_orig(ref1,1)) (x_orig(ref1)+10000*v_orig(ref1,1))],[(y_orig(ref1)-10000*v_orig(ref1,2)) (y_orig(ref1)+10000*v_orig(ref1,2))],[(z_orig(ref1)-10000*v_orig(ref1,3)) (z_orig(ref1)+10000*v_orig(ref1,3))],'Color','r','LineWidth',1);
end
xlim([aoi(1)-10,aoi(2)+10]);
ylim([aoi(3)-10,aoi(4)+10]);
zlim([aoi(5)-10,aoi(6)+10]);
view(0,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  User Defined Bisector Reference Point  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Prompt User to Define Point for X & Z References
% fprintf('Using intuition, "Click" on figure(1) approximate cast-off center of origin region.\n'); %User Instructions
% fprintf('"Click" ENTER when complete.\n\n'); %User Instructions
% time1 = toc;
% [x_ref_xz, z_ref_xz] = getpts; %Select on figure(1) Intuitive Cast-off Center of Origin Region
% time2 = toc;
% clc %Remove User Instructions 
% Ref_xz = [x_ref_xz, z_ref_xz]; %Create Reference Point
% 
% %Plot Room Dimensions, Stains, and Stain Trajectories
% figure(2)
% hold on 
% grid on
% axis equal
% % axis(room_size)
% h1_front_yz = plot([aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','c','LineWidth',5); %Plot Front Surface Dimensions
% h1_downward_yz = plot([aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(5) aoi(5) aoi(5) aoi(5) aoi(5)],'Color','g','LineWidth',4); %Plot Downward Surface Dimensions
% h1_back_yz = plot([aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','b','LineWidth',3); %Plot Back Surface Dimensions
% h1_upward_yz = plot([aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(6) aoi(6) aoi(6) aoi(6) aoi(6)],'Color','y','LineWidth',2); %Plot Downward Surface Dimensions
% h2_yz = plot(Ys,Zs, '.','MarkerSize',max_room_size*0.1,'Color','r','LineWidth',2); %Plot XYZ Stains
% title({'Cast-off Center & Radius of Origin',''})
% xlabel(['Y-Axis (cm)'])
% ylabel(['Z-Axis (cm)'])
% 
% for ref2 = 1:numstains
%     h3_yz(ref2) = plot([(Ys(ref2)-1000*v(ref2,2)) (Ys(ref2)+1000*v(ref2,2))],[(Zs(ref2)-1000*v(ref2,3)) (Zs(ref2)+1000*v(ref2,3))]);
% end
% h_yz = [h1_front_yz h1_downward_yz h1_back_yz h1_upward_yz h2_yz h3_yz(1)]; % h9(sample2) h10
% legend(h_yz, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Trajectories')
% 
% %Prompt User to Define Point for Y Reference
% fprintf('Using intuition, "Click" on figure(2) approximate cast-off center of origin region. If cast-off circle is not apparent choose approximate stain trajectory clustered centroid.\n'); %User Instructions
% fprintf('"Click" ENTER when complete.\n\n'); %User Instructions
% time3 = toc;
% [y_ref_yz, z_ref_yz] = getpts; %Select on figure(2) Intuitive Cast-off Center of Origin Region
% time4 = toc;
% clc %Remove User Instructions and Clear Figure
% 
% Ref_yz = [y_ref_yz, z_ref_yz]; %Create Reference Point
% 
% hold off
% 
% Ref = [x_ref_xz y_ref_yz z_ref_xz]; %Compiled User Defined Refernce Point
% x_ref = Ref(1); %X-coordinate of User Defined Refernce Point
% y_ref = Ref(2); %Y-coordinate of User Defined Refernce Point
% z_ref = Ref(3); %Z-coordinate of User Defined Refernce Point

% stains = 1:length(Xs); %Stain Element Reference
% nums = downsample(stains,dwnsamp); %Downsample Stains as User Defined
% num = nums(1:sampsize); %Select Cluster Group by Sample Size

clstr = linspace(1,numstains,numstains); %Possible Clusters
clstr_comb = combvec(combvec(clstr,clstr),clstr); %All Possible Combinations and Permutations of Clusters
clstr_comb(:,any(diff(sort(clstr_comb,1),[],1)==0,1))=[]; %Remove Repeated Combinations
[~,iclstr,~] = unique(sort(clstr_comb',2),'rows'); %Determine Repeated Permutations
clstr_comb = clstr_comb(:,sort(iclstr))'; %Remove Repeated Permutations

clstr_num = size(clstr_comb,1); %Number of Clusters

if user_clstr == 1;
    clstr_num = size(stain_cluster,1); %Determine Number of Clusters from User Defined Clusters for Analysis
end

anti_v = [zeros(size(v,1),1) zeros(size(v,1),1) 1*ones(size(v,1),1)]; %Antigravity Vector
an = cross(v,anti_v,2); %Normal Vector of Plane defined by Velocity Vector and Antigravity Vector
ann = sqrt(sum(an.^2,2)); %Normalization of Normal Vector
for ain = 1:numstains;
    an(ain,:) = an(ain,:)/ann(ain);
end
au = cross([zeros(numstains,1) -1*ones(numstains,1) zeros(numstains,1)],an,2); %Determine Line of Intersection between XZ-plane and Plane Best Fitting Stain Velocity Vectors
anu = sqrt(sum(au.^2,2)); %Normal of Best Plane 
au1 = au./[anu anu anu];
%     au1 = [zeros(numstains,1) zeros(numstains,1) -1*ones(numstains,1)];
aphi1 = -acos(dot([zeros(numstains,1) -1*ones(numstains,1) zeros(numstains,1)],an,2)); %Planar Angle to Rotate Plane about Line of Intersection to XZ-plane

for akn = 1:numstains; % 'aphi1' Corrections
    if aphi1(akn) < (-1*0.5*pi)
        aphi1(akn) = aphi1(akn)+pi;
    elseif aphi1(akn) > (0.5*pi)
        aphi1(akn) = aphi1(akn)-pi;
    end
end

Vn = v; %Save Original Velocity Vectors
xs = Xs; %Save Original X-coordinates
ys = Ys; %Save Original Y-coordinates
zs = Zs; %Save Original Z-coordinates
alpha_orig = alpha; %Save Original Alpha Angles
beta_orig = beta; %Save Original Beta Angles
    
    for rain = 1:numstains;
        ar1(:,:,rain) = [((au1(rain,1)^2)+((au1(rain,2)^2)+(au1(rain,3)^2))*cos(aphi1(rain))) (au1(rain,1)*au1(rain,2)*(1-cos(aphi1(rain)))-au1(rain,3)*sin(aphi1(rain))) (au1(rain,1)*au1(rain,3)*(1-cos(aphi1(rain)))+au1(rain,2)*sin(aphi1(rain))); (au1(rain,1)*au1(rain,2)*(1-cos(aphi1(rain)))+au1(rain,3)*sin(aphi1(rain))) ((au1(rain,2)^2)+((au1(rain,1)^2)+(au1(rain,3)^2))*cos(aphi1(rain))) (au1(rain,2)*au1(rain,3)*(1-cos(aphi1(rain)))-au1(1)*sin(aphi1(rain))); (au1(rain,1)*au1(rain,3)*(1-cos(aphi1(rain)))-au1(rain,2)*sin(aphi1(rain))) (au1(rain,2)*au1(rain,3)*(1-cos(aphi1(rain)))+au1(1)*sin(aphi1(rain))) ((au1(rain,3)^2)+((au1(rain,1)^2)+(au1(rain,2)^2))*cos(aphi1(rain)))]; %3D Rotation Matrix
        vproj(rain,:) = (ar1(:,:,rain)*v(rain,:)')'; %Projected Velocity Vector
        PG0(rain,:) = [Xs(rain) Zs(rain)]; %Initial X and Z-coordinates
        PG1(rain,:) = [Xs(rain) (Zs(rain)+10)]; %X and Z-coordinates of Antigravity Line
        PG2(rain,:) = [(Xs(rain)+vproj(rain,1)) (Zs(rain)+vproj(rain,3))]; %X and Z-coordinates plus Projected Vector
        pgn1(rain,:) = (PG2(rain,:)-PG0(rain,:))/norm(PG2(rain,:)-PG0(rain,:));
        pgn2(rain,:) = (PG1(rain,:)-PG0(rain,:))/norm(PG1(rain,:)-PG0(rain,:));
        alpha_pg(rain,:) = abs(atan2(norm(det([pgn2(rain,:);pgn1(rain,:)])),dot(pgn1(rain,:),pgn2(rain,:)))); %Projected Global Alpha angle relative to Gravity
        dot_tv(rain,:) = dot(e_t(rain,:),vproj(rain,:),2); %Directionality Check of Projected Vectors

        %Projected Global Alpha Angle Surface Corrections
        if face(rain) == 1;
            if dot_tv(rain)<=0;
                alpha_gn(rain,:) = 2*pi-alpha_pg(rain);
            else;
                alpha_gn(rain,:) = 2*pi-alpha_pg(rain);
            end;
        elseif face(rain) == 2;
            if dot_tv(rain)<=0;
                alpha_gn(rain,:) = 2*pi-alpha_pg(rain);
            else;
                alpha_gn(rain,:) = alpha_pg(rain);
            end;
            alpha(rain,:) = (0.5*pi)-alpha_pg(rain);
        elseif face(rain) == 3;
            if dot_tv(rain)<=0;
                alpha_gn(rain,:) = alpha_pg(rain);
            else;
                alpha_gn(rain,:) = alpha_pg(rain);
            end
        elseif face(rain) == 4;
            if dot_tv(rain)<=0;
                alpha_gn(rain,:) = 2*pi-alpha_pg(rain);
            else;
                alpha_gn(rain,:) = alpha_pg(rain);
            end
        end
    end

time5 = toc();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Clustering by Global Alpha Impact Angle  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if alpha_clstr == 1;
    for ac = 1:length(clstr_alpha);
        Alpha_gn = combvec(alpha_gn',alpha_gn')'; %All Possible Combinations of Two Global Alpha Angles
        ALPHA_GN = combvec(Alpha_gn',Alpha_gn(:,1)')'; %All Possible Combinations of Three Global Alpha Angles
        pgn_ind = combvec(1:numstains,1:numstains)'; %All Possible Combinations of Two Stain Indices
        PGN_IND = combvec(pgn_ind',pgn_ind(:,1)')'; %All Possible Combinations of Three Stain Indices

        ALPHA_GN(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        PGN_IND(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        [~,idk,~] = unique(sort(PGN_IND,2),'rows'); %Determine Repeated Permutations
        ALPHA_GN = ALPHA_GN(sort(idk),:); %Remove Repeated Permutations
        PGN_IND = PGN_IND(sort(idk),:); %Remove Repeated Permutations

        for io = 1:size(PGN_IND,1);
            if abs(ALPHA_GN(io,1)-ALPHA_GN(io,2))>=(clstr_alpha(ac,:)-alpha_std) && ...
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,2))<=(clstr_alpha(ac,:)+alpha_std) && ... 
                    abs(ALPHA_GN(io,2)-ALPHA_GN(io,3))>=(clstr_alpha(ac,:)-alpha_std) && ... 
                    abs(ALPHA_GN(io,2)-ALPHA_GN(io,3))<=(clstr_alpha(ac,:)+alpha_std) || ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,3))>=(clstr_alpha(ac,:)-alpha_std) && ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,3))<=(clstr_alpha(ac,:)+alpha_std) && ... 
                    abs(ALPHA_GN(io,2)-ALPHA_GN(io,3))>=(clstr_alpha(ac,:)-alpha_std) && ... 
                    abs(ALPHA_GN(io,2)-ALPHA_GN(io,3))<=(clstr_alpha(ac,:)+alpha_std) || ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,2))>=(clstr_alpha(ac,:)-alpha_std) && ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,2))<=(clstr_alpha(ac,:)+alpha_std) && ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,3))>=(clstr_alpha(ac,:)-alpha_std) && ... 
                    abs(ALPHA_GN(io,1)-ALPHA_GN(io,3))<=(clstr_alpha(ac,:)+alpha_std);
                comb_mat(io,:,ac) = PGN_IND(io,:);
            else
                comb_mat(io,:,ac) = [NaN, NaN, NaN]; %Remove Errored Clusters
            end
        end

        clstr_num(ac,:) = sum((sum(~isnan(comb_mat(:,:,ac)),2)==3),1); %Total Number of Clusters for this Clustering Method
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Clustering by Distance between Stains  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dist_clstr == 1 && dwn_samp_stains ~= 1;
    for ad = 1:length(clstr_dist);
        r_dist = [xs ys zs]; %Stain Coordinates
        R_dist = combvec(r_dist',r_dist')'; %All Possible Combinations and Permutations of Two Sets of Stain Coordinates
        R_DIST = combvec(R_dist',R_dist(:,1:3)')'; %All Possible Combinations and Permutations of Three Sets of Stain Coordinates
        DIST_12 = sqrt(sum((R_DIST(:,1:3)-R_DIST(:,4:6)).^2,2)); %Distance between First Two Stains of All Possible Combinations and Permutations of Stains
        DIST_23 = sqrt(sum((R_DIST(:,4:6)-R_DIST(:,7:9)).^2,2)); %Distance between Last Two Stains of All Possible Combinations and Permutations of Stains
        DIST_13 = sqrt(sum((R_DIST(:,1:3)-R_DIST(:,7:9)).^2,2)); %Distance between Middle Two Stains of All Possible Combinations and Permutations of Stains
        pgn_ind = combvec(1:numstains,1:numstains)'; %All Possible Combinations and Permutations of Two Sets of Stain Indices
        PGN_IND = combvec(pgn_ind',pgn_ind(:,1)')'; %All Possible Combinations and Permutations of Three Sets of Stain Indices

        DIST_12(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        DIST_23(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        DIST_13(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        PGN_IND(any(diff(sort(PGN_IND,2),[],2)==0,2),:) = []; %Remove Repeated Combinations
        [~,idl,~] = unique(sort(PGN_IND,2),'rows'); %Determine Repeated Permutations
        DIST_12 = DIST_12(sort(idl),:); %Remove Repeated Permutations
        DIST_23 = DIST_23(sort(idl),:); %Remove Repeated Permutations
        DIST_13 = DIST_13(sort(idl),:); %Remove Repeated Permutations
        PGN_IND = PGN_IND(sort(idl),:); %Remove Repeated Permutations

        for in = 1:size(PGN_IND,1);
            if DIST_12(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_12(in) <= (clstr_dist(ad,:)+dist_std) && ...
                    DIST_23(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_23(in) <= (clstr_dist(ad,:)+dist_std) || ...
                    DIST_13(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_13(in) <= (clstr_dist(ad,:)+dist_std) && ...
                    DIST_23(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_23(in) <= (clstr_dist(ad,:)+dist_std) || ...
                    DIST_12(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_12(in) <= (clstr_dist(ad,:)+dist_std) && ...
                    DIST_13(in) >= (clstr_dist(ad,:)-dist_std) && ...
                    DIST_13(in) <= (clstr_dist(ad,:)+dist_std);
                comb_mat(in,:,ad) = PGN_IND(in,:);
            else
                comb_mat(in,:,ad) = [NaN, NaN, NaN]; %Remove Errored Clusters
            end
        end

        clstr_num(ad,:) = sum((sum(~isnan(comb_mat(:,:,ad)),2)==3),1); %Total Number of Clusters for this Clustering Method
    end
elseif dist_clstr == 1 && dwn_samp_stains == 1;
    for it = 1:length(dwnsamp);
        clstr_num(it,:) = ((numstains-((sampsize-1)*dwnsamp(it)))/(A(it)-overlap(it))); %Number of Clusters
        for iv = 1:clstr_num(it);
            dwn(iv,:) = ((1:dwnsamp(it):sampsize*dwnsamp(it)))+(A(it)-overlap(it))*(iv-1); %Select Downsampling Criteria
            r_dist(iv,:) = [xs(dwn(iv,1)) ys(dwn(iv,1)) zs(dwn(iv,1)) xs(dwn(iv,2)) ys(dwn(iv,2)) zs(dwn(iv,2)) xs(dwn(iv,3)) ys(dwn(iv,3)) zs(dwn(iv,3))]; %Select Stian Coordinates
            dist_12(iv,:) = sqrt(sum((r_dist(iv,1:3)-r_dist(iv,4:6)).^2,2)); %Determine Distance between First Two Clustered Stains
            dist_23(iv,:) = sqrt(sum((r_dist(iv,4:6)-r_dist(iv,7:9)).^2,2)); %Determine Distance between Last Two Clustered Stains
            dist_13(iv,:) = sqrt(sum((r_dist(iv,1:3)-r_dist(iv,7:9)).^2,2)); %Determine Distance between Middle Two Clustered Stains
        end
        comb_mat = dwn; %Total Number of Clusters for this Clustering Method
        R_DIST = r_dist; %Total Distance between Stains
        DIST_12 = dist_12; %Total Distance between First Two Clustered Stains
        DIST_23 = dist_23; %Total Distance between Last Two Clustered Stains
        DIST_13 = dist_13; %Total Distance between Middle Two Clustered Stains
        Comb_mat = dwn; %Total Number of Clusters for this Clustering Method
        clearvars dwn r_dist dist_12 dist_23 dist_13 %Remove Variables to avoid Overwrite Errors
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Cluster Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if alpha_clstr == 1 || dist_clstr == 1;
    comb_num = size(comb_mat,3); %Determine Total Number of Stains for These Clustering Methods
else
    comb_num = 1;  %Determine Total Number of Stains for All Possible Clustered Stain Combinations
end

results = cell(comb_num,9); %Save Results

for iq = 1:comb_num;
    if alpha_clstr == 1 || dist_clstr == 1 && dwn_samp_stains ~= 1;
        comb_mat_1 = comb_mat(:,1,iq); %Determine Total Number of Stains for the First Two Stains in these Clustering Methods
        comb_mat_2 = comb_mat(:,2,iq); %Determine Total Number of Stains for the Last Two Stains in these Clustering Methods
        comb_mat_3 = comb_mat(:,3,iq); %Determine Total Number of Stains for the Middle Two Stains in these Clustering Methods
        Comb_mat{iq,:} = [comb_mat_1(~isnan(comb_mat_1(:,:))) comb_mat_2(~isnan(comb_mat_2(:,:))) comb_mat_3(~isnan(comb_mat_3(:,:)))]; %Determine Total Number of Stains for These Clustering Methods

        if isempty(Comb_mat{iq}) ==1; %Skip Cluster if Empty
            continue
        end
    end
    
    isocubes = ones(numel(cu_cx),1); %Pre-allocate PDF Distribution Values
    
    for ip = 1:clstr_num(iq); %Number of Clusters
        if user_clstr == 1;
            B(ip,:) = stain_cluster(ip,:); %Cluster Stain Indices
        else
            if alpha_clstr == 1 || dist_clstr == 1;
                B = Comb_mat; %Cluster Stain Indices
            elseif dwn_samp_stains == 1 && opti_space ~= 1;
                B(ip,:) = ((1:dwnsamp:sampsize*dwnsamp))+(A-overlap)*(ip-1); %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i>1
            else
                B = clstr_comb; %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i>1
            end
        end
        V = v(B(ip,:),:); %Clustered Velocity Vector
        XS = Xs(B(ip,:),:); %Clustered X-coordinate of Stain Location (with Gravity and Drag Contributions if Selected)
        YS = Ys(B(ip,:),:); %Clustered Y-coordinate of Stain Location (with Gravity and Drag Contributions if Selected)
        ZS = Zs(B(ip,:),:); %Clustered Z-coordinate of Stain Location (with Gravity and Drag Contributions if Selected)
        x_s = xs(B(ip,:),:); %Clustered X-coordinate of Stain Location (without Gravity and Drag Contributions)
        y_s = ys(B(ip,:),:); %Clustered Y-coordinate of Stain Location (without Gravity and Drag Contributions)
        z_s = zs(B(ip,:),:); %Clustered Z-coordinate of Stain Location (without Gravity and Drag Contributions)
        Alpha_p = alpha_p(B(ip,:),:); %Clustered Projected Alpha
        Alpha = alpha(B(ip,:),:); %Clustered Alpha (with Gravity and Drag Contributions if Selected)
        Alpha_pg = alpha_pg(B(ip,:),:); %Clustered Global Projected Alpha
        Alpha_orig = alpha_orig(B(ip,:),:); %Clustered Alpha (without Gravity and Drag Contributions)
        Gamma = gamma(B(ip,:),:); %Clustered Gamma
        Minor = minor(B(ip,:),:); %Clustered Minor Axes
        Face = face(B(ip,:),:); %Clustered Stain Face of Impact
        E_n = e_n(B(ip,:),:); %Clustered Normal Tangential Vectors
        E_t = e_t(B(ip,:),:); %Clustered Tangential Vectors
        E_nxt = e_nxt(B(ip,:),:); %Clustered Normal x Tangential Vectors
        [Weight,Sn] = Castoff_Reconstruction_FUNC(Face,V,aoi,XS,YS,ZS,Alpha_p,Alpha,Alpha_pg,Alpha_orig,Gamma,Minor,Ref,InOutTrajectory,InRoom,max_room_size,min_room_size,Lx,Ly,Lz,Nx,Ny,Nz,res,xmin,ymin,zmin,stdev,cu_cx,cu_cy,cu_cz,ip,isocubes,Spread_Fact_cu,Spread_Fact_theta,Spread_Fact_upsilon,dalpha_range,dgamma_range,datamat,clstr_num,iq,comb_num); %Function Determining Castoff Reconstruction
        figure(3);
        hold on;
        plot3(xs(B(ip,:)),ys(B(ip,:)),zs(B(ip,:)),'.','MarkerSize',max_room_size*0.1,'Color','g','LineWidth',2);
        quiver3(xs(B(ip,:)),ys(B(ip,:)),zs(B(ip,:)),v(B(ip,:),1),v(B(ip,:),2),v(B(ip,:),3),'Color','g');

        V_s(ip,:,:) = V;
        X_S(ip,:) = XS';
        Y_S(ip,:) = YS';
        Z_S(ip,:) = ZS';
        X_s(ip,:) = x_s';
        Y_s(ip,:) = y_s';
        Z_s(ip,:) = z_s';
        AlphA(ip,:) = Alpha';
        Alpha_PG(ip,:) = Alpha_pg';
        Alpha_ORIG(ip,:) = Alpha_orig';
        GammA(ip,:) = Gamma';
        FacE(ip,:) = Face';

        if ip == 1
            if any(Weight) == 0
                isocubes = ones(size(Weight)); %Skip PDF Distribution of All Zeros
            else
                isocubes = isocubes.*Weight; %Add PDF Distribution for first Iteration
            end
        else
            if any(Weight) == 0
                isocubes = isocubes; %Skip PDF Distribution of All Zeros
            else
                isocubes = isocubes.*Weight; %Multiply PDF Distribution for Remaining Iterations
            end
        end
    
    if ip == 1;
        NUM_clstr = clstr_num; %Count Initial PDF Distribution
    else
        if any(Weight == 0);
            NUM_clstr = NUM_clstr - 1; %Skip PDF Distributions of All Zeros
        end
    end
    
    end

    figure(3);
    hold off;
    xlim([aoi(1)-10,aoi(2)+10]);
    ylim([aoi(3)-10,aoi(4)+10]);
    zlim([aoi(5)-10,aoi(6)+10]);
    
%     if (dwn_samp_stains == 1 && (dist_clstr == 1 || alpha_clstr == 1)) == 1;
%         results{iq,1} = [DIST_12{:,:,iq} DIST_23{:,:,iq} DIST_13{:,:,iq}]; %Save Results
%     elseif dist_clstr == 1 && dist_clstr == 1;
%         results{iq,1} = clstr_dist(iq); %Save Results
%     elseif alpha_clstr == 1;
%         results{iq,1} = clstr_alpha(iq)*180/pi; %Save Results
%     else
%         results{iq,1} = iq; %Save Results
%     end
%     results{iq,2} = B; %Save Results
%     results{iq,3} = [X_S Y_S Z_S]; %Save Results
%     results{iq,4} = [X_s Y_s Z_s]; %Save Results
%     results{iq,5} = AlphA; %Save Results
%     results{iq,6} = Alpha_PG; %Save Results
%     results{iq,7} = Alpha_ORIG; %Save Results
%     results{iq,8} = GammA; %Save Results
%     results{iq,9} = FacE; %Save Results
   
    clearvars V_s X_S Y_S Z_S X_s Y_s Z_s AlphA Alpha_PG Alpha_ORIG GammA FacE cgrade %Remove Variables to avoid Overwriting Errors
    
    isocubes = reshape(isocubes,[Nz,Nx,Ny]); %Reshape Array to 3D Cube Orientation
    maxcu = max(max(max(isocubes,[],1),[],2),[],3); %Determine Maximum Distributed Product Value
    isonorm = isocubes/maxcu; %Normalize Product Distribution
    isoscale = (isocubes(:)-min(isocubes(:)))./(max(isocubes(:))-min(isocubes(:))); %Scale Product Distributed Values between 0 and 1
    isoscale = reshape(isoscale,[Nz,Nx,Ny]); %Reshape Array to 3D Cube Orientation
    isoscale(isoscale==0)=realmin('double');

    [f1,v1] = isosurface(Xcu,Ycu,Zcu,isonorm,percentiles(1)^NUM_clstr); %Save Distribution Vertices of First Percentile
    [f2,v2] = isosurface(Xcu,Ycu,Zcu,isonorm,percentiles(2)^NUM_clstr); %Save Distribution Vertices of Second Percentile
    [f3,v3] = isosurface(Xcu,Ycu,Zcu,isonorm,percentiles(3)^NUM_clstr); %Save Distribution Vertices of Third Percentile

time6 = toc;

%Plot Resultant Scaled Product Distribution Percentiles
fig_num = 3+iq;
figure(fig_num);
hold on;
grid on;
p1_front = plot3([aoi(1) aoi(1) aoi(1) aoi(1) aoi(1)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','c','LineWidth',5); %Plot Front Surface Dimensions
p1_downward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(5) aoi(5) aoi(5) aoi(5) aoi(5)],'Color','g','LineWidth',4); %Plot Downward Surface Dimensions
p1_back = plot3([aoi(2) aoi(2) aoi(2) aoi(2) aoi(2)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','b','LineWidth',3); %Plot Back Surface Dimensions
p1_upward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(6) aoi(6) aoi(6) aoi(6) aoi(6)],'Color','y','LineWidth',2); %Plot Downward Surface Dimensions
p9 = plot3(x_orig,y_orig,z_orig,'.','Color','c','MarkerSize',10); %Plot Original XYZ Stains
p2 = plot3(Xs,Ys,Zs,'.','Color','r','MarkerSize',15); %Plot XYZ Stains
chi = linspace(0, 2*pi, 25); %Angle Vector
% ************************************************
% % ***** Use when Cast-off Motion is known *****
% x_actual = actual_x + actual_r.*cos(chi); %Resultant X-Coordinate of Cast-off Circle
% y_actual = actual_y*ones(size(x_actual)); %Resultant Y-Coordinate of Cast-off Circle
% z_actual = actual_z + actual_r.*sin(chi); %Resultant Z-Coordinate of Cast-off Circle
% ************************************************
for ref1 = 1:numstains;
    p10(ref1) = plot3([(x_orig(ref1)-10000*v_orig(ref1,1)) (x_orig(ref1)+10000*v_orig(ref1,1))],[(y_orig(ref1)-10000*v_orig(ref1,2)) (y_orig(ref1)+10000*v_orig(ref1,2))],[(z_orig(ref1)-10000*v_orig(ref1,3)) (z_orig(ref1)+10000*v_orig(ref1,3))],'Color','c','LineWidth',1);
end;
for ref1 = 1:numstains;
    p8(ref1) = plot3([(Xs(ref1)-10000*v(ref1,1)) (Xs(ref1)+10000*v(ref1,1))],[(Ys(ref1)-10000*v(ref1,2)) (Ys(ref1)+10000*v(ref1,2))],[(Zs(ref1)-10000*v(ref1,3)) (Zs(ref1)+10000*v(ref1,3))],'Color','r','LineWidth',1);
end;
% ************************************************
% % ***** Use when Cast-off Motion is known *****
% p3 = plot3(x_actual,y_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
% p4 = plot3(actual_x,actual_y,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
% ************************************************
colorcu = {[1,0,0];[0,1,0];[0,0,1];[0,1,1]};
transcu = [1.0 0.5 0.2 0.2 0.1 0.075 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];

if isempty(nonzeros(isonorm-ones(size(isonorm)))) == 1;
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
else
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)^NUM_clstr)),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   p6 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)^NUM_clstr)),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
   p7 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(3)^NUM_clstr)),'FaceColor',cell2mat(colorcu(3)),'EdgeAlpha',transcu(4),'FaceAlpha',transcu(3));
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p6 p7 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(percentiles(1)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(percentiles(2)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(percentiles(3)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
end

title({'Castoff Reconstruction'});
xlabel(['X-Axis (cm)']);
ylabel(['Y-Axis (cm)']);
zlabel(['Z-Axis (cm)']);
view(-30,30);
set(gca,'FontSize',20);
axis equal;
% set(gcf, 'Position', get(0, 'Screensize')); %Make Figure Full-screen
hold off;
xlim([aoi(1)-10,aoi(2)+10]);
ylim([aoi(3)-10,aoi(4)+10]);
zlim([aoi(5)-10,aoi(6)+10]);
% fig_3 = figure(3);
% close(fig_3);
end

save(datamat); %Save Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Output Results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(regexprep(datamat,'.mat','_OUTPUT.txt'),'w'); %Open OUTPUT
fprintf(fileID, '%s\r\n', 'Cast-off Reconstruction Results:'); %Display OUTPUT Variables
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Faces(x)','High Percentile Faces(y)','High Percentile Faces(z)'); %Display OUTPUT Variables
for outi = 1:length(f1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f1(outi,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Vertices(x)','High Percentile Vertices(y)','High Percentile Vertices(z)'); %Display OUTPUT Variables
for outj = 1:length(v1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v1(outj,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','Medium Percentile Faces(x)','Medium Percentile Faces(y)','Medium Percentile Faces(z)'); %Display OUTPUT Variables
for outk = 1:length(f2)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f2(outk,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','Medium Percentile Vertices(x)','Medium Percentile Vertices(y)','Medium Percentile Vertices(z)'); %Display OUTPUT Variables
for outl = 1:length(v2)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v2(outl,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Faces(x)','High Percentile Faces(y)','High Percentile Faces(z)'); %Display OUTPUT Variables
for outm = 1:length(f3)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f3(outm,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Vertices(x)','High Percentile Vertices(y)','High Percentile Vertices(z)'); %Display OUTPUT Variables
for outn = 1:length(v3)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v3(outn,:)'); %Output Results to OUTPUT.txt
end
fclose(fileID);

fprintf('Script ran at %s\n', datestr(now,'HH:MM:SS.FFF')); %End Run Designator with Time
tot_time = toc();
fprintf('Total Elapsed Run Time: %s\n', tot_time);
fprintf('Total Elapsed Cluster Analysis Time: %s\n', (time6-time5));
min_W_dist = min(min(min(isocubes,[],1),[],2),[],3);
fprintf('Minimum Weight Distribution: %4.2d', min_W_dist);
fprintf(' ,  Note: if "Minimum Weight Distribution" is zero, consider including fewer clusters to prevent distributing zeros %s\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%  Send Email Notification of Completion  %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%  (use for long iterative compilations)  %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datatxt = regexprep(datamat,'^','COMPLETED Castoff Reconstruction for ');
% setpref('Internet','SMTP_Server','mailhub.iastate.edu');
% setpref('Internet','E_mail','scottres@iastate.edu');
% sendmail('scottres@iastate.edu',datatxt);
