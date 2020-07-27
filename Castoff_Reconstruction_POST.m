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
% % % 2. Set desired DRIVER (lines 150-269) and MAIN (lines 134-158) User Inputs along with Ad-hoc User defined Variables. (FARO and Hemospat Drivers are provided for trialing code)
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
% % %  - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (lines 157-159)
% % %  - 'lngth','minor' major and minor axis lengths of stains in mm (lines 164-165)
% % %  - 'alpha','gamma' stain impact and directional angles (lines 166-168)
% % %  - 'room', room length (along x-dimension), width (along y-dimension), and height (along z-direction) (x,y,z) in centimeters (line 199)
% % %  - 'xmin','xmax','ymin','ymax','zmin','zmax' minimum and maximum coordinate of possible region of cast-off origin) (lines 193-199)
% % %  - 'res' Spatial resolution of reconstruction (Length of Discretized Uniform Regions of Space Dimensions) (1-15cm is the recommended range) (15cm took ~30 seconds, 10cm took ~1 minute, and 7.5cm took ~1 hour in a large room with several hundred stains and default specifications) (line 203)
% % %  - 'n_b','n_u','n_f','n_d','t_b','t_u','t_f','t_d' surface normal (n) and tangential (t) unit vectors (four surfaces by default are assumed to be perpendicular) (surfaces include: back,upward,front,downward) *** Does not incorporate side surfaces *** (lines 223-230)
% % % - 'inclsurf_[1-4]' Choose '1' to INCLUDE Surface #[1-4] Stains; Choose '0' to EXCLUDE Surface #[1-4]
% % %  - 'InOutTrajectory' *****NOT RECOMMENDED TO USE***** Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
% % %  - 'alpha60more' Choose '1' to Remove impact angles Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees (line 254)
% % %  - 'alphaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Impact Angles (alpha); Choose '0' to NOT Weight Center & Radius Results by Impact Angles (alpha)
% % %  - 'gammaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Glancing Angles (gamma); Choose '0' to NOT Weight Center & Radius Results by Glancing Angles (gamma)
% % %  - 'zetaweight' *****ALWAYS USE == 1***** Choose '1' to Weight Center & Radius Results by Bisector Intersection Angles (zeta); Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta)
% % %  - 'bsctintweight' *****NOT RECOMMENDED TO USE***** Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta) (line 259)
% % %  - 'InRoom' *****NOT RECOMMENDED TO USE***** Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories
% % % 
% % % User Inputs: (MAIN)
% % %  - 'percentiles' Choose Resultant Product Distribution Percentiles for Plotting Likelihood Regions (this can be change when replotting the final figure (Figure 4) in the Post-processor (line 133)
% % % 
% % % Ad-hoc User Defined Variables:
% % %  - 'stdev' Approximate Standard Deviation of Stain Width Measurement in mm (0.1 is the default value) (DRIVER line 201)
% % %  - 'Spread_Fact_cu' Spreading Factor Uncertainty in Distance between a given Discretized Region in Space and Arc (10cm-30cm is recommended range, 2*res (cm) is the default) (DRIVER line 237)
% % %  - 'Spread_Fact_theta' Spreading Factor Uncertainty in In-plane Angle (Theta) in radians (pi/180 radians (1 degree) is the default) (DRIVER line 238)
% % %  - 'Spread_Fact_upsilon' Spreading Factor Uncertainty in Off-plane Angle (Upsilon) in radians (pi/18-pi/6 radians (10-30 degrees) is recommended pi/9 radians (20 degrees) is the default (DRIVER line 239)
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

%Note to User:
%This is a Post Processor for 'Castoff_Reconstruction_MAIN.m'

clear,clc,close all

load('Ink_Trial_DRIVER.mat');

percentiles = [1e-1, 1e-4, 1e-7];

figure(4);
hold on;
grid on;
p1_front = plot3([aoi(1) aoi(1) aoi(1) aoi(1) aoi(1)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','c','LineWidth',5); %Plot Front Surface Dimensions
p1_downward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(5) aoi(5) aoi(5) aoi(5) aoi(5)],'Color','g','LineWidth',4); %Plot Downward Surface Dimensions
p1_back = plot3([aoi(2) aoi(2) aoi(2) aoi(2) aoi(2)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','b','LineWidth',3); %Plot Back Surface Dimensions
p1_upward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(6) aoi(6) aoi(6) aoi(6) aoi(6)],'Color','y','LineWidth',2); %Plot Downward Surface Dimensions
p2 = plot3(Xs,Ys,Zs, '.','MarkerSize',max_room_size*0.1,'Color','r','LineWidth',2); %Plot XYZ Stains
chi = linspace(0, 2*pi, 25); %Angle Vector
x_actual = actual_x + actual_r.*cos(chi); %Resultant X-Coordinate of Cast-off Circle
y_actual = actual_y*ones(size(x_actual)); %Resultant Y-Coordinate of Cast-off Circle
z_actual = actual_z + actual_r.*sin(chi); %Resultant Z-Coordinate of Cast-off Circle
for ref1 = 1:numstains;
    p8(ref1) = plot3([(Xs(ref1)-10000*v(ref1,1)) (Xs(ref1)+10000*v(ref1,1))],[(Ys(ref1)-10000*v(ref1,2)) (Ys(ref1)+10000*v(ref1,2))],[(Zs(ref1)-10000*v(ref1,3)) (Zs(ref1)+10000*v(ref1,3))],'Color','r','LineWidth',1);
end;
p3 = plot3(x_actual,y_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
p4 = plot3(actual_x,actual_y,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
colorcu = {[1,0,0];[0,1,0];[0,0,1];[0,1,1]};
transcu = [1.0 0.5 0.2 0.2 0.1 0.075 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];
p5 = patch(isosurface(Xcu,Ycu,Zcu,isoscale,percentiles(1)),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
p6 = patch(isosurface(Xcu,Ycu,Zcu,isoscale,percentiles(2)),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
p7 = patch(isosurface(Xcu,Ycu,Zcu,isoscale,percentiles(3)),'FaceColor',cell2mat(colorcu(3)),'EdgeAlpha',transcu(4),'FaceAlpha',transcu(3));
p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p3 p4 p5 p6 p7];
legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', 'Actual Castoff Circle Location', 'Actual Castoff Center Location', 'high Percentile Castoff Reconstruction', 'mid Percentile Castoff Reconstruction', 'low Percentile Castoff Reconstruction');
title({'Castoff Reconstruction'});
xlabel(['X-Axis (cm)']);
ylabel(['Y-Axis (cm)']);
zlabel(['Z-Axis (cm)']);
view(-30,30);
axis square;
xlim([min_room_size-100,max_room_size+100]);
ylim([min_room_size-100,max_room_size+100]);
zlim([min_room_size-100,max_room_size+100]);
set(gcf, 'Position', get(0, 'Screensize')); %Make Figure Full-screen
hold off;
