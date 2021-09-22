% % % %%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
% % % Reconstructs stains from cast-off event to reproduce the motion of cast-off.
% % % Last Updated 02/01/2021
% % % 
% % % Licenses:
% % % All licenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the repository Licenses folder, this was not intentional by the author.
% % % 
% % % To install and run code, follow instructions in README.txt
% % %
% % % User Inputs: (INPUT & DRIVER)
% % %  - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (INPUT)
% % %  - 'lngth','minor' major and minor axis lengths of stains in mm (INPUT)
% % %  - 'alpha','gamma' stain impact and directional angles (INPUT)
% % %  - 'room', room length (along x-dimension), width (along y-dimension), and height (along z-direction) (x,y,z) in centimeters (INPUT)
% % %  - 'xmin','xmax','ymin','ymax','zmin','zmax' minimum and maximum coordinate of possible region of cast-off origin) (DRIVER)
% % %  - 'res' Spatial resolution of reconstruction (Length of Discretized Uniform Regions of Space Dimensions) (1-15cm is the recommended range) (15cm took ~30 seconds, 10cm took ~1 minute, and 7.5cm took ~1 hour in a large room with several hundred stains and default specifications) (INPUT)
% % %  - 'n_b','n_u','n_f','n_d','n_l','n_r','t_b','t_u','t_f','t_d','t_l','t_r' surface normal (n) and tangential (t) unit vectors (six surfaces by default are assumed to be perpendicular) (surfaces include: back,upward,front,downward,left,right) (DRIVER)
% % %  - 'inclsurf_[1-6]' Choose '1' to INCLUDE Surface #[1-6] Stains; Choose '0' to EXCLUDE Surface #[1-6]
% % %  - 'InOutTrajectory' *****NOT RECOMMENDED TO USE***** Choose '1' to Select Stain Trajectories ONLY Directed into Room Dimensions; Choose '0' to Select Trajectories Directed BOTH In and Out of Room Dimensions
% % %  - 'alpha30less' Choose '1' to Remove impact angles Alpha Values Less All Alpha Values Greater Than 30 degrees
% % %  - 'alpha60more' Choose '1' to Remove impact angles Alpha Values Greater Than 60 degrees; Choose '0' to Keep All Alpha Values Greater Than 60 degrees (DRIVER)
% % %  - 'bsctintweight' *****NOT RECOMMENDED TO USE***** Choose '1' to Weight Center & Radius Results by Distance between Stain Trajectory Intersection & Bisector Intersection, d; Choose '0' to NOT Weight Center & Radius Results by Bisector Intersection Angles (zeta) (DRIVER)
% % %  - 'InRoom' *****NOT RECOMMENDED TO USE***** Choose '1' to Remove Bisector Trajectories Outside of Room Dimensions; Choose '0' to Keep All Bisector Trajectories (DRIVER)
% % % 
% % % User Inputs: (MAIN)
% % %  - 'percentiles' Choose Resultant Product Distribution Percentiles for Plotting Likelihood Regions (this can be change when replotting the final figure (Figure 4) in the Post-processor (MAIN)
% % % 
% % % Ad-hoc User Defined Variables:
% % %  - 'stdev' Approximate Standard Deviation of Stain Width Measurement in mm (0.1 is the default value) (INPUT)
% % %  - 'Spread_Fact_cu' Spreading Factor Uncertainty in Distance between a given Discretized Region in Space and Arc (10cm-30cm is recommended range, 2*res (cm) is the default) (MAIN)
% % %  - 'Spread_Fact_theta' Spreading Factor Uncertainty in In-plane Angle (Theta) in radians (pi/180 radians (1 degree) is the default) (MAIN)
% % %  - 'Spread_Fact_upsilon' Spreading Factor Uncertainty in Off-plane Angle (Upsilon) in radians (pi/18-pi/6 radians (10-30 degrees) is recommended pi/9 radians (20 degrees) is the default (MAIN)
% % %  
% % % Recommended Clustering Methods:
% % %  - 'stain_cluster' Enter user defined stains indices in an nx3 matrix for n clusters to cluster specific stain combinations (MAIN)
% % %  - 'dwn_samp_stains' Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option (MAIN)
% % %  - 'dwnsamp' Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains (MAIN)
% % %  - 'sampsize' Select Stain Cluster Sample Size by Integer Factor Greater than Three (3) (MAIN)
% % %  - 'overlap' Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired; Set Equal to 'dwnsamp*(sampsize-1)' for maximum overlap (MAIN)
% % %  - 'opti_space' Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1). (MAIN)

%Note to User:
%This is a Post Processor for 'Castoff_Reconstruction_MAIN.m'

clear,clc,close all

load('INK_Trial_INPUT.mat');

number_traces = 1; %Number of Cast-off Swings
regions = [0.95, 0.75,0.60]
percentiles = [regions(1)^(number_traces*NUM_clstr), regions(2)^(number_traces*NUM_clstr), regions(3)^(number_traces*NUM_clstr)]

fig_num = 3+iq;
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
colorcu = {[1,0,0];[0,1,0];[0,0,1];[0,1,1]};
transcu = [1.0 0.5 0.2 0.2 0.1 0.075 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];

if isempty(nonzeros(isonorm-ones(size(isonorm)))) == 1 || ~(any(isonorm(:) < (percentiles(1))));
   p = [p1_front p1_downward p1_back p1_upward p2 p3 p4 p8(1) p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Actual Castoff Circle Location', 'Actual Castoff Center Location', 'Stain Straight-line Trajectories', 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
elseif ~(any(isonorm(:) < (percentiles(2))));
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   volumes_cc = abs([meshVolume(v1,[],f1)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)))*res^3
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(regions(1)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
elseif ~(any(isonorm(:) < (percentiles(3))));
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   p6 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)),'noshare'),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
   volumes_cc = abs([meshVolume(v1,[],f1) meshVolume(v2,[],f2)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)))*res^3
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p6 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(regions(1)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(regions(2)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
else   
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   p6 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)),'noshare'),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
   p7 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(3)),'noshare'),'FaceColor',cell2mat(colorcu(3)),'EdgeAlpha',transcu(4),'FaceAlpha',transcu(3));
   volumes_cc = abs([meshVolume(v1,f1) meshVolume(v2,f2) meshVolume(v3,f3)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)))*res^3
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p6 p7 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(regions(1)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(regions(2)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(regions(3)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
end

title({'Castoff Reconstruction'});
xlabel(['X-Axis (cm)']);
ylabel(['Y-Axis (cm)']);
zlabel(['Z-Axis (cm)']);
view(S_n);
set(gca,'FontSize',20);
axis equal;
hold off;
xlim([aoi(1)-10,aoi(2)+10]);
ylim([aoi(3)-10,aoi(4)+10]);
zlim([aoi(5)-10,aoi(6)+10]);
