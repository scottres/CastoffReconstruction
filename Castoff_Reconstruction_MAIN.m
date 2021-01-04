% % % %%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
% % % Reconstructs stains from cast-off event to reproduce the motion of cast-off.
% % % Last Updated 11/30/2020
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
% % %  - 'n_b','n_u','n_f','n_d','t_b','t_u','t_f','t_d' surface normal (n) and tangential (t) unit vectors (four surfaces by default are assumed to be perpendicular) (surfaces include: back,upward,front,downward) *** Does not incorporate side surfaces *** (DRIVER)
% % %  - 'inclsurf_[1-4]' Choose '1' to INCLUDE Surface #[1-4] Stains; Choose '0' to EXCLUDE Surface #[1-4]
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  User Defined Values  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = 'INK_Trial_INPUT_DRIVER.mat'; %Producing figure(3) and figure(4+) which simulate Figures 5 and 7 of 'McCleary et al FSI 2021', respectively OR data = 'FARO_Trial_10_INPUT_DRIVER.mat'; %Producing figure(3) and figure(4+) which simulate Figures 5 and 6 of 'McCleary et al FSI 2021', respectively

load(data); %'Castoff_Reconstruction_DRIVER_SWINEBLOOD.mat' OR 'Castoff_Reconstruction_DRIVER_LISCIO_TRIAL_2_UPDATED.mat'
sig_n = 3; %Minimum PDF Value Standard Deviations
percentiles = [0.90 0.75 0.60]; %Choose Resultant Product Distribution Percentiles for Plotting (this can be change when replotting the final figure (Figure 4)
SF_cu_range = [0,400]; %Distance Spreading Factor Allowable Range, in cm (zero to 5cm or 30 cm give reasonable results)
dalpha_range = [0,0.5*pi]; %Alpha Impact Angle Range
dgamma_range = [0,2*pi]; %Gamma Directional Angle Range
res_range = [0,15]; %Spatial Region Resolution Allowable Range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datamat = regexprep(data,'_DRIVER','','ignorecase'); %Change .mat name for Saving Results
dlmwrite(regexprep(datamat,'.mat','_OUTPUT.txt'),[]); %Create Output .txt-file

dwn_samp_stains = 1; %Set Equal to '1' to Cluster by Downsampling, Set Equal to '0' for another Cluster Option
dwnsamp = 1; %[1:15]'; %Select Stain Cluster Sample Rate by Integer Factor; Enter '1' if Clustering Adjacent Stains
sampsize = 3; %round(0.45*numstains); %Select Stain Cluster Sample Size by Integer Factor Greater than Three (3)
overlap = dwnsamp*(sampsize-1); %sampsize-1; %Select Cluster Sampling Overlap by Integer Factor; Enter '0' if No Overlap is Desired
A = 2*dwnsamp+(sampsize-2); %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i=1

opti_space = 1; %Set Equal to '1' to Analyze Spatter with Equally Spaced Stains (by Distance between Stains), Set Equal to '0' to Not Apply Equal Spacing. Use with Downsampling (dwn_samp_stains = 1).

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
if res >= res_range(1) && res < res_range(2)
  fprintf(fileID, '%s\r\n', strcat('res =','',num2str(res))); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "res" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('res =','',num2str(res),' (ideal "res" Range =','',num2str(res_range))); %Inputs are within Predefined Range
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
if opti_space == opti_default
  fprintf(fileID, '%s\r\n%s\r\n', strcat('opti_space =','',num2str(opti_space)),''); %Inputs are within Predefined Range
else
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "opti_space" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n%s\r\n', strcat('opti_space =','',num2str(opti_space),' (ideal "opti_space" value =','',num2str(opti_default),')'),''); %Inputs are within Predefined Range
end

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
v_n = [e_n(:,1).*v_n e_n(:,2).*v_n e_n(:,3).*v_n]; %Normal Component of Velocity Vector
v_nxt = [e_nxt(:,1).*v_nxt e_nxt(:,2).*v_nxt e_nxt(:,3).*v_nxt]; %Normal x Tangential Component of Velocity Vector
v_t = [e_t(:,1).*v_t e_t(:,2).*v_t e_t(:,3).*v_t]; %Tangential Component of Velocity Vector
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Automatically Determine Bisector Reference Point  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine Reference Point within Region of Cast-off Origin
xzi = [Xs Ys Zs];%[x y z]+max_room_size*v; %Initial Trajectory Endpoint
xzf = [Xs Ys Zs]-max_room_size*v; %Final Trajectory Endpoint
xrefi = combvec2(xzi(:,1)',xzi(:,1)')'; %Compiling All Possible Trajectory Final X-Coordinate Combinations and Purmutations
zrefi = combvec2(xzi(:,3)',xzi(:,3)')'; %Compiling All Possible Trajectory Final Z-Coordinate Combinations and Purmutations
xreff = combvec2(xzf(:,1)',xzf(:,1)')'; %Compiling All Possible Trajectory Final X-Coordinate Combinations and Purmutations
zreff = combvec2(xzf(:,3)',xzf(:,3)')'; %Compiling All Possible Trajectory Final Z-Coordinate Combinations and Purmutations
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
chi = linspace(0, 2*pi, 500); %Angle Vector for Plotting Known Cast-off Motion

% ************************************************
% % ***** Use when Cast-off Motion is known *****
%  x_fin = (center(1) + radius.*cos(chi))'; %Resultant X-Coordinate of Cast-off Circle
%  y_fin = (center(2)*ones(size(x_fin)));
%  z_fin = (center(3) + radius.*sin(chi))'; %Resultant Z-Coordinate of Cast-off Circle

 x_actual = actual_x + actual_r.*cos(chi); %Resultant X-Coordinate of Cast-off Circle
 y_actual = actual_y*ones(size(x_actual)); %Resultant Y-Coordinate of Cast-off Circle
 z_actual = actual_z + actual_r.*sin(chi); %Resultant Z-Coordinate of Cast-off Circle

 h15 = plot(x_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
 h16 = plot(actual_x,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
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
% ************************************************
% % ***** Use when Cast-off Motion is known *****
%  x_fin = (center(1) + radius.*cos(chi))'; %Resultant X-Coordinate of Cast-off Circle
%  y_fin = (center(2)*ones(size(x_fin)));
%  z_fin = (center(3) + radius.*sin(chi))'; %Resultant Z-Coordinate of Cast-off Circle

 x_actual = actual_x + actual_r.*cos(chi); %Resultant X-Coordinate of Cast-off Circle
 y_actual = actual_y*ones(size(x_actual)); %Resultant Y-Coordinate of Cast-off Circle
 z_actual = actual_z + actual_r.*sin(chi); %Resultant Z-Coordinate of Cast-off Circle

 h15 = plot3(x_actual,y_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
 h16 = plot3(actual_x,actual_y,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
% ************************************************
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
clstr_comb = combvec2(combvec2(clstr,clstr),clstr); %All Possible Combinations and Permutations of Clusters
clstr_comb(:,any(diff(sort(clstr_comb,1),[],1)==0,1))=[]; %Remove Repeated Combinations
[~,iclstr,~] = unique(sort(clstr_comb',2),'rows'); %Determine Repeated Permutations
clstr_comb = clstr_comb(:,sort(iclstr))'; %Remove Repeated Permutations

clstr_num = size(clstr_comb,1); %Number of Clusters

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

for akn = 1:numstains; %Planar Angle Corrections
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

if dwn_samp_stains == 1;
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

comb_num = 1;  %Determine Total Number of Stains for All Possible Clustered Stain Combinations

for iq = 1:comb_num;
        
    isocubes = ones(numel(cu_cx),1); %Pre-allocate PDF Distribution Values
    
    for ip = 1:clstr_num(iq); %Number of Clusters
        if dwn_samp_stains == 1;
            B = Comb_mat; %Cluster Stain Indices
        elseif dwn_samp_stains == 1 && opti_space ~= 1;
            B(ip,:) = ((1:dwnsamp:sampsize*dwnsamp))+(A-overlap)*(ip-1); %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i>1
        else
            B = clstr_comb; %Equation Defining Relationship between dwnsamp, sampsize, and overlap for i>1
        end
        V = v(B(ip,:),:); %Clustered Velocity Vector
        XS = Xs(B(ip,:),:); %Clustered X-coordinate of Stain Location
        YS = Ys(B(ip,:),:); %Clustered Y-coordinate of Stain Location
        ZS = Zs(B(ip,:),:); %Clustered Z-coordinate of Stain Location
        x_s = xs(B(ip,:),:); %Clustered X-coordinate of Stain Location
        y_s = ys(B(ip,:),:); %Clustered Y-coordinate of Stain Location
        z_s = zs(B(ip,:),:); %Clustered Z-coordinate of Stain Location
        Alpha_p = alpha_p(B(ip,:),:); %Clustered Projected Alpha
        Alpha = alpha(B(ip,:),:); %Clustered Alpha
        Alpha_pg = alpha_pg(B(ip,:),:); %Clustered Global Projected Alpha
        Alpha_orig = alpha_orig(B(ip,:),:); %Clustered Alpha
        Gamma = gamma(B(ip,:),:); %Clustered Gamma
        Minor = minor(B(ip,:),:); %Clustered Minor Axes
        Face = face(B(ip,:),:); %Clustered Stain Face of Impact
        E_n = e_n(B(ip,:),:); %Clustered Normal Tangential Vectors
        E_t = e_t(B(ip,:),:); %Clustered Tangential Vectors
        E_nxt = e_nxt(B(ip,:),:); %Clustered Normal x Tangential Vectors
        [Weight,Sn] = Castoff_Reconstruction_FUNC(Face,V,aoi,XS,YS,ZS,Alpha_p,Alpha,Alpha_pg,Alpha_orig,Gamma,Minor,Ref,InOutTrajectory,InRoom,max_room_size,min_room_size,Lx,Ly,Lz,Nx,Ny,Nz,res,xmin,ymin,zmin,stdev,cu_cx,cu_cy,cu_cz,ip,isocubes,DeltaGamma(B(ip,:),:),dalpha_range,dgamma_range,SF_cu_range,datamat,clstr_num,iq,comb_num,sig_n,fileID); %Function Determining Castoff Reconstruction
        if ishandle(3)
          figure(3);
          hold on;
          plot3(xs(B(ip,:)),ys(B(ip,:)),zs(B(ip,:)),'.','MarkerSize',max_room_size*0.1,'Color','g','LineWidth',2);
          quiver3(xs(B(ip,:)),ys(B(ip,:)),zs(B(ip,:)),v(B(ip,:),1),v(B(ip,:),2),v(B(ip,:),3),'Color','g');
          xlim([aoi(1)-10,aoi(2)+10]);
          ylim([aoi(3)-10,aoi(4)+10]);
          zlim([aoi(5)-10,aoi(6)+10]);
        end

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

        scaling=log10(realmax('double'))-max(0,log10(max(isocubes)))-max(0,log10(max(Weight)))-5; %Rescale joint PDF results to avoid realmin('double') and realmax('double')
        isocubes=isocubes*10^scaling;

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
       
    percent_calculation = strcat(num2str(round(ip/clstr_num(iq)*100)),'%') %Display Completion Precentage
    Current_min = min(min(min(isocubes,[],1),[],2),[],3); %Display Current Minimum
    if ip == 1;
        NUM_clstr = clstr_num; %Count Initial PDF Distribution
    else
        if any(Weight == 0);
            NUM_clstr = NUM_clstr - 1; %Skip PDF Distributions of All Zeros
        end
    end
end
  
   hold off;
   
    clearvars V_s X_S Y_S Z_S X_s Y_s Z_s AlphA Alpha_PG Alpha_ORIG GammA FacE cgrade %Remove Variables to avoid Overwriting Errors
    
    isocubes = reshape(isocubes,[Nz,Nx,Ny]); %Reshape Array to 3D Cube Orientation
    maxcu = max(max(max(isocubes,[],1),[],2),[],3); %Determine Maximum Distributed Product Value
    isonorm = isocubes/maxcu; %Normalize Product Distribution
    isoscale = (isocubes(:)-min(isocubes(:)))./(max(isocubes(:))-min(isocubes(:))); %Scale Product Distributed Values between 0 and 1
    isoscale = reshape(isoscale,[Nz,Nx,Ny]); %Reshape Array to 3D Cube Orientation
    isoscale(isoscale==0)=realmin('double');

    [f1,v1] = isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)^NUM_clstr),'noshare'); %Save Distribution Vertices of First Percentile
    [f2,v2] = isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)^NUM_clstr),'noshare'); %Save Distribution Vertices of Second Percentile
    [f3,v3] = isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(3)^NUM_clstr),'noshare'); %Save Distribution Vertices of Third Percentile

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
 p3 = plot3(x_actual,y_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
 p4 = plot3(actual_x,actual_y,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
% ************************************************
colorcu = {[1,0,0];[0,1,0];[0,0,1];[0,1,1]};
transcu = [1.0 0.5 0.2 0.2 0.1 0.075 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];

if isempty(nonzeros(isonorm-ones(size(isonorm)))) == 1 || ~(any(isonorm(:) < (percentiles(1)^NUM_clstr)));
   p = [p1_front p1_downward p1_back p1_upward p2 p3 p4 p8(1) p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Actual Castoff Circle Location', 'Actual Castoff Center Location', 'Stain Straight-line Trajectories', 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
elseif ~(any(isonorm(:) < (percentiles(2)^NUM_clstr)));
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   volumes_cc = abs([meshVolume(v1,[],f1)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)^NUM_clstr))*res^3
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(percentiles(1)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
elseif ~(any(isonorm(:) < (percentiles(3)^NUM_clstr)));
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   p6 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
   volumes_cc = abs([meshVolume(v1,[],f1) meshVolume(v2,[],f2)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)^NUM_clstr))*res^3
   p = [p1_front p1_downward p1_back p1_upward p2 p8(1) p5 p6 p9 p10(1)]; % p3 p4
   legend(p, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Straight-line Trajectories', strcat(num2str(percentiles(1)*100), 'th Percentile Castoff Reconstruction'), strcat(num2str(percentiles(2)*100), 'th Percentile Castoff Reconstruction'), 'Spatter Stains not Included in Analysis', 'Stain Straight-line Trajectories not Included in Analysis', 'Location', 'northeastoutside'); % 'Actual Castoff Circle Location', 'Actual Castoff Center Location',
else   
   p5 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(1)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(1)),'EdgeAlpha',transcu(1),'FaceAlpha',transcu(1));
   p6 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(2)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(2)),'EdgeAlpha',transcu(3),'FaceAlpha',transcu(2));
   p7 = patch(isosurface(Xcu,Ycu,Zcu,isonorm,(percentiles(3)^NUM_clstr),'noshare'),'FaceColor',cell2mat(colorcu(3)),'EdgeAlpha',transcu(4),'FaceAlpha',transcu(3));
   volumes_cc = abs([meshVolume(v1,f1) meshVolume(v2,f2) meshVolume(v3,f3)]) %Vol_Region_1 = (sum(isonorm(:) > percentiles(1)^NUM_clstr))*res^3
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
% fileID = fopen(regexprep(datamat,'.mat','_OUTPUT.txt'),'w'); %Open OUTPUT
fprintf(fileID, '%s\r\n', 'Cast-off Reconstruction Results:'); %Display OUTPUT Variables
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Faces(x)','High Percentile Faces(y)','High Percentile Faces(z)'); %Display OUTPUT Variables
for outi = 1:size(f1,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f1(outi,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Vertices(x)','High Percentile Vertices(y)','High Percentile Vertices(z)'); %Display OUTPUT Variables
for outj = 1:size(v1,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v1(outj,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','Medium Percentile Faces(x)','Medium Percentile Faces(y)','Medium Percentile Faces(z)'); %Display OUTPUT Variables
for outk = 1:size(f2,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f2(outk,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','Medium Percentile Vertices(x)','Medium Percentile Vertices(y)','Medium Percentile Vertices(z)'); %Display OUTPUT Variables
for outl = 1:size(v2,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v2(outl,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Faces(x)','High Percentile Faces(y)','High Percentile Faces(z)'); %Display OUTPUT Variables
for outm = 1:size(f3,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',f3(outm,:)'); %Output Results to OUTPUT.txt
end
fprintf(fileID, '%12s %24s %36s\r\n','High Percentile Vertices(x)','High Percentile Vertices(y)','High Percentile Vertices(z)'); %Display OUTPUT Variables
for outn = 1:size(v3,1)
  fprintf(fileID, '%12.3f %24.3f %36.3f\r\n',v3(outn,:)'); %Output Results to OUTPUT.txt
end
fclose(fileID);

fprintf('Script ran at %s\n', datestr(now,'HH:MM:SS.FFF')); %End Run Designator with Time
tot_time = toc();
fprintf('Total Elapsed Run Time: %s\n', strcat(num2str(tot_time), ' seconds'));
fprintf('Total Elapsed Cluster Analysis Time: %s\n', strcat(num2str(time6-time5), ' seconds'));
min_W_dist = min(min(min(isocubes,[],1),[],2),[],3);
fprintf('Minimum Weight Distribution: %s', num2str(min_W_dist));
fprintf(',  Note: if "Minimum Weight Distribution" is zero, consider increasing spreading factors, increasing "sig_n", or including fewer clusters to prevent distributing zeros %s\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%  Send Email Notification of Completion  %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%  (use for long iterative compilations)  %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datatxt = regexprep(datamat,'^','COMPLETED Castoff Reconstruction for ');
% setpref('Internet','SMTP_Server','your_email_server_here');
% setpref('Internet','E_mail','your_email_address_here');
% sendmail('your_email_address_here',datatxt);
