%%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
Reconstructs stains from cast-off event to reproduce the motion of cast-off.

%Primary Developer: Scott McCleary; Supervisor: Daniel Attinger
%Emails: scott.thomas.mccleary@gmail.com | daniel.attinger@gmail.com
%Phone: (515) 975-5544
%Cast-off Spatter Reconstruction

%D.A. and S.MC acknowledge financial support from the Center for Statistics and Applications in Forensic Evidence (CSAFE), the National Institute of Justice (NIJ), and Iowa State University, as well as pro bono technical and consulting services from Struo LLC, a scientific consulting company based in Ames IA.  D.A. thanks his Iowa State University colleague Xuan Hien Nguyen for pointing the geometry proposition to him.  E.L. would like to thank Independent Forensic Services (IFS) for their assistance and use of their facilities and resources in collecting the non-circular cast-off patterns.
%Primary Developer: Scott McCleary; Supervisor: Daniel Attinger
%Emails: scott.thomas.mccleary@gmail.com | daniel.attinger@gmail.com
%Phone: (515) 975-5544
%Spatter Stains to Cast-off  Spatter Reconstruction
%Work supported by Center for Statistics and Applications in Forensic Evidence
%Department of Mechanical Engineering - Attinger Lab
%Iowa State University
%04/28/2020

Cast-off Reconstruction has only been tested on Matlab 2019a and GNU Octave 5.2.0. It is recommended that Matlab 2019a, GNU Octave 5.2.0 or newer versions are used to run the included files. Users run the included files at their own risk. Especially, D.A. and S.MC assume no responsibility regarding the results, their interpretation and the consequences thereof. 

Required Repository Files to run the code:
  - Spatter Measurement Data, e.g. 'Ink_Trial.csv', 'Swineblood_Trial.csv'
- 'Castoff_Reconstruction_DRIVER.m'
 	 - 'DRIVER.csv' produces 'DRIVER.mat' required for  'Castoff_Reconstruction_MAIN.m'
MAIN
 - 'Castoff_Reconstruction_MAIN.m'
 - 'Castoff_Reconstruction_FUNC.m'
 - 'lineSegmentIntersect.m'
 - 'point_to_line.m'
 - 'gauss_distribution.m'
 - 'CircleFitByPratt.m'
 - Desired Data, e.g. 'Ink_Trial.csv', 'Swineblood_Trial.csv'
 - 'Castoff_Reconstruction_POST.m' for Post-processing
 - 'meshVolume.m' for Post-processor
 - 'inpolyhedron.m' for Post-processor
 - 'triangulateFaces.m' for Post-processor

Liscenses:
All liscenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the included documentationrepository ‘Licenses’ folder, this was not intentional by the author.

Required Installation for GNU Octave Compatibility: (Suggested Installation Order)
 - Updated GNU Octave (written and tested with 5.2.0) (https://www.gnu.org/software/octave/download.html)
	- All packages can be updated to the latest version by running:
  		>> pkg update
 - Updated Java (https://www.java.com/en/download/win10.jsp)
 - Psychtoolbox Toolbox Version 3 (PTB-3) (http://psychtoolbox.org/)
	Psychtoolbox-3 Requirements (http://psychtoolbox.org/download.html#download-problems)
	 - Subversion 1.7.x command-line client (e.g. Sliksvn) (https://sliksvn.com/download/)
	 - 64-Bit GStreamer-1.16.0 MSVC runtime or later versions (https://gstreamer.freedesktop.org/data/pkg/windows/1.16.0/gstreamer-1.0-msvc-x86_64-1.16.0.msi)
 - 'combvec.m'

Instructions to Run:
1. Save all included files to the same directory.
2. Set desired  DRIVER (lines 110-229) and MAIN (lines 94-117) User Inputs along with Ad-hoc User defined Variables.
3. Choose desired clustering method ('dwn_samp_stains' and 'opti_space' is  the default 'Best' Clustering Method).
4. 'Run' (F5) the DRIVER ('Castoff_Reconstruction_DRIVER.m').
5. 'Run' (F5) the MAIN ('Castoff_Reconstruction_MAIN.m').
6. MAIN will output the 'Total Elapsed Cluster Analysis Time:', in seconds (s), with the total program run time excluding any variable time from user input. See ‘User Inputs: (DRIVER)’ ‘res’ for estimated runtimes.
7. MAIN will save the Resultant Variables as a .mat file denoted by ‘Castoff_Reconstruction.mat’ and Output the Cast-off Reconstruction Results as Figures (4+) (figure number is dependent on clustering method).

User Inputs: (DRIVER)
 - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (lines 117-119)
 - 'lngth','minor' major and minor axis lengths of stains in mm (lines 124-125)
 - 'alpha','gamma' stain pitching impact and glancing impactdirectional  angles (lines 126-128)
 - 'room', room dimensions length (along x-dimension), width (along y-dimension), and height (along z-direction) (x,y,z) in centimeters (line 159)
 - 'xmin','xmax','ymin','ymax','zmin','zmax' minimum and maximum surface coordinates for area of interestcoordinate of (plausiblepossible region of cast-off origin) (lines 153-159)
 - 'res' Spatial Rresolution of Heat Mapreconstruction (Length of Discretized Uniform Regions of Space Cube Dimensions) (1-15cm is the recommended range) (15cm took ~30 seconds, 10cm took ~1 minute, and 7.5cm took ~1 hour in a large room with several hundred stains and default specifications) (line 163)
 - 'n_b','n_u','n_f','n_d','t_b','t_u','t_f','t_d' surface normal (n) and tangential (t) unit vectors (four surfaces by default surfaces are plumbed andassumed to be perpendicular) (surfaces include: back,upward,front,downward) *** Does not incorporate side surfaces *** (lines 183-190)
- 'inclsurf_[1-4]' Choose '1' to INCLUDE Surface #[1-4] Stains; Choose '0' to EXCLUDE Surface #[1-4]
 - 'InOutTrajectory' *****NOT RECOMMENDED TO USE***** Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
 - 'alpha60more' Choose '1' to Remove impact angles Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees (line 214)
 - 'alphaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Impact Angles (alpha); Choose '0' to NOT Weight Center & Radius Results by Impact Angles (alpha)
 - 'gammaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Glancing Angles (gamma); Choose '0' to NOT Weight Center & Radius Results by Glancing Angles (gamma)
 - 'zetaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Bisector Intersection Angles (zeta); Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
 - 'bsctintweight' *****NOT RECOMMENDED TO USE***** Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta) (line 219)
 - 'InRoom' *****NOT RECOMMENDED TO USE***** Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories

User Inputs: (MAIN)
 - 'percentiles' Choose Resultant Product Distribution Percentiles for Plotting Likelihood Regions (this can be change when replotting the final figure (Figure 4) in the Post-processor (line 94)

Ad-hoc User Defined Variables:
 - 'stdev' Approximated Standard Deviation of Stain Width Measurement in mm (0.1 is the default value) (DRIVER line 161)
 - 'Spread_Fact_cu' Spreading Factor Uncertainty in Distance between a given Discretized Region in Space Cube and Arc in centimeters (10cm-30cm is recommended range, 2*res (cm) is the default) (DRIVER line 197)
 - 'Spread_Fact_theta' Spreading Factor Uncertainty in In-plane Angle (Theta) in radians (pi/180 radians (1 degree) is the default) (DRIVER line 198)
 - 'Spread_Fact_upsilon' Spreading Factor Uncertainty in Off-plane Angle (Upsilon) in radians (pi/18-pi/6 radians (10-30 degrees) is recommended pi/9 radians (20 degrees) is the default (DRIVER line 199)
 
Recommended Clustering Methods:
 - 'user_clstr' Set Equal to '1' to Run One Specific Cluster of Three Stains, Set Equal to '0' for another Cluster Option
 - 'stain_cluster' Enter user defined stains indices in an nx3 matrix for n clusters to cluster specific stain combinations
 - 'dwn_samp_stains' Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option (MAIN line 102)
 - 'dwnsamp' Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains (MAIN line 103)
 - 'sampsize' Select Stain Cluster Sample Size by Integer Factor Greater than Three (3) (MAIN line 104)
 - 'overlap' Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired; Set Equal to 'dwnsamp*(sampsize-1)' for maximum overlap (MAIN line 105)
 - 'alpha_clstr' Set Equal to '1' to Cluster by Half Global Alpha Impact Angle, Set Equal to '0' for another Cluster Option
 - 'clstr_alpha' Select Difference in Half Global Alpha Impact Angle for Clustering in radians
 - 'alpha_std' Select Standard Deviation of Difference in Half Global Alpha Impact Angle for Clustering in radians
 - 'dist_clstr' Set Equal to '1' to Cluster by Distance between Stains, Set Equal to '0' for another Cluster Option
 - 'clstr_dist' Select Distance between Stains for Clustering in centimeters
 - 'dist_std' Select Standard Deviation of Distance between Stains for Clustering in centimeters
 - 'opti_space' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1). (MAIN line 116)
 - 'opti_angle' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Angle between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).

Instructions to Run:
1. Save all included files to the same directory.
2. Set desired DRIVER and MAIN User Inputs along with Ad-hoc User defined Variables.
3. Choose desired clustering method ('dwn_samp_stains' and 'opti_space' is the default 'Best' Clustering Method).
4. 'Run' (F5) the DRIVER ('Castoff_Reconstruction_DRIVER.m').
5. 'Run' (F5) the MAIN ('Castoff_Reconstruction_MAIN.m').
6. MAIN will output the 'Total Elapsed Cluster Analysis Time:', in seconds (s), with the total program run time excluding any variable time from user input. 
7. MAIN will output the Cast-off Reconstruction Results in as Figures (4+) (figure number is dependantdependent on clustering method).

