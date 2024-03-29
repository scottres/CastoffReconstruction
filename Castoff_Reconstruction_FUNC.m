function [Psi_tot,Sn] = Castoff_Reconstruction_FUNC(face,v,aoi,Xs,Ys,Zs,alpha_p,alpha,alpha_pg,alpha_orig,gamma,minor,Ref,InOutTrajectory,InRoom,max_room_size,min_room_size,Lx,Ly,Lz,Nx,Ny,Nz,res,xmin,ymin,zmin,stdev,cu_cx,cu_cy,cu_cz,ip,isocubes,deltagamma,dalpha_range,dgamma_range,SF_cu_range,datamat,clstr_num,iq,comb_num,sig_n,fileID,S_n); %Function Determining Castoff Reconstruction
% % % %%%%% MATLAB/Octave Cast-off Reconstruction %%%%%
% % % Reconstructs stains from cast-off event to reproduce the motion of cast-off.
% % %
% Scott McCleary
% Email: scott.thomas.mccleary@gmail.com | daniel.attinger@gmail.com
% Phone: (515) 975-5544
% Spatter Stains to Cast-off Reconstruction
% Center for Statistics and Applications in Forensic Evidence
% Department of Mechanical Engineering - Attinger Lab
% Iowa State University
% % %
% % % Last Updated 09/23/2021
% % % 
% % % Required Repository Files to run the code:
% % %  - Spatter Measurement Data, e.g. 'Ink_Trial_INPUT.csv', 'Swineblood_Trial_INPUT.csv'
% % %  - 'Castoff_Reconstruction_DRIVER.m'
% % %  - 'DRIVER.csv' produces 'DRIVER.mat' required for 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_MAIN.m'
% % %  - 'Castoff_Reconstruction_FUNC.m'
% % %  - 'lineSegmentIntersect.m'
% % %  - 'meshVolume.m'
% % %  - 'point_to_line.m'
% % %  - 'gauss_distribution.m'
% % %  - 'CircleFitByPratt.m'
% % %  - 'Castoff_Reconstruction_POST.m'
% % %  - 'combvec2.m'
% % %  - 'inpolyhedron.m' 
% % %  - 'triangulateFaces.m' 
% % %  - 'linecirc.m'
% % %  - 'generate_input.m'
% % %    - 'linecirc.m'
% % %    - 'plane_line_intersect.m'
% % %    - 'rotate_3D.m'
% % %    - 'triangulateFaces.m'
% % %   
% % % Licenses:
% % % All licenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the repository Licenses folder, this was not intentional by the author.
% % % 
% % % To install and run code, follow instructions in README.txt
% % %
% % % User Inputs: (INPUT & DRIVER)
% % %  - 'x','y','z' stain coordinate locations (in centimeters) relative to user defined origin (INPUT)
% % %  - 'lngth','minor' major and minor axis lengths of stains in mm (INPUT)
% % %  - 'gamma' stain impact and directional angles (INPUT)
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

numstains = size(Xs,1); %Recalculate Number of Stains being Clustered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Determine Best Plane  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

limit1s = [aoi(1) aoi(3) aoi(5); aoi(2) aoi(3) aoi(5); aoi(2) aoi(3) aoi(6); aoi(1) aoi(3) aoi(6); aoi(1) aoi(3) aoi(5)]; %Keep Track of Left Room Dimensions throughout Translations & Rotations
limit2s = [aoi(1) aoi(4) aoi(5); aoi(2) aoi(4) aoi(5); aoi(2) aoi(4) aoi(6); aoi(1) aoi(4) aoi(6); aoi(1) aoi(4) aoi(5)]; %Keep Track of Right Room Dimensions throughout Translations & Rotations

%Find Normal Vector of Plane that Best Fits Given Elements in V
Xn0 = S_n; %Initial Guess/Starting Point for 'fminsearch' Iteration
options = optimset('MaxIter',1e10,'Algorithm','levenberg-marquardt','Display','off','MaxFunEvals',1e5); %Option to View 'fsolve' Iteration
funn = @(Xn)sum((v*[Xn(1),Xn(2),Xn(3)]').^2); %Function to Minimize Dot Product of All Stain Velocities and Choosen Normal Vector
Xn = fsolve(funn,Xn0,options); %Minimization of Function, Add ",options" to Procedure to View 'fminsearch' Iteration
norm_n = sqrt(sum([Xn(1),Xn(2),Xn(3)].^2,2)); %Normal of Minimization Function Result
Sn = [Xn./norm_n]; %Normalize Minimization Function Result
if dot(Sn,Xn0)<0
    Sn = -Sn;
end

%Find Point Closest to the Given Stains
Xp0 = [0.5*(aoi(2)-aoi(1)),0.5*(aoi(4)-aoi(3)),0.5*(aoi(6)-aoi(5))]; %Initial Guess/Starting Point for 'fminsearch' Iteration 
funp = @(Xp)sum(sqrt(sum(([Xs Ys Zs]-[Xp(1).*ones(numstains,1) Xp(2).*ones(numstains,1) Xp(3).*ones(numstains,1)]).^2))); %Function to Minimize Distance between Stain Locations and Choosen Point
Xp = fminsearch(funp,Xp0); %Minimization of Function, Add ",options" to Procedure to View 'fminsearch' Iteration

%Project Stain Velocity Vectors to Best Plane
subXYZp = (([Xs,Ys,Zs]-[Xp(1).*ones(numstains,1) Xp(2).*ones(numstains,1) Xp(3).*ones(numstains,1)])*Sn'); %Subsidary Step to Determine Project Stains to Best Plane
XYZp = [Xs Ys Zs]-[subXYZp.*Sn(1) subXYZp.*Sn(2) subXYZp.*Sn(3)]; %Project Stains to Best Plane
Vp = v-[((v*Sn')./(sqrt(sum(Sn.^2,2)).^2)).*Sn(1) ((v*Sn')./(sqrt(sum(Sn.^2,2)).^2)).*Sn(2) ((v*Sn')./(sqrt(sum(Sn.^2,2)).^2)).*Sn(3)]; %Project Stain Velocities to Best Plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  2D to 3D Translation and Rotation  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Translate Stains to Origin
d = Sn*Xp'; %Dot Product of Normal Vector and Point Closest to the Given Stains
pt2pln = [((d*Sn(1))./sqrt(sum(Sn.^2))) ((d*Sn(2))./sqrt(sum(Sn.^2))) ((d*Sn(3))./sum(sqrt(Sn.^2)))]; %Shortest Distance between Origin and Plane
XYZt = [(XYZp(:,1)-Xp(1).*ones(numstains,1)) (XYZp(:,2)-Xp(2).*ones(numstains,1)) (XYZp(:,3)-Xp(3).*ones(numstains,1))]; %Apply Translation of Projected Stains to Origin by Point Closest to Given Stains

%Rotate Stains and Stain Vectors to XZ-plane
u = cross([0,-1,0],Sn); %Determine Line of Intersection between XZ-plane and Plane Best Fitting Stain Velocity Vectors
norm_u = sqrt(sum([u(1),u(2),u(3)].^2,2)); %Normal of Best Plane
U1 = [u./norm_u]; %Normalized Best Plane
phi1 = -acos(dot([0,-1,0],Sn)); %Planar Angle to Rotate Plane about Line of Intersection to XZ-plane

R1 = [((U1(1)^2)+((U1(2)^2)+(U1(3)^2))*cos(phi1)) (U1(1)*U1(2)*(1-cos(phi1))-U1(3)*sin(phi1)) (U1(1)*U1(3)*(1-cos(phi1))+U1(2)*sin(phi1)); (U1(1)*U1(2)*(1-cos(phi1))+U1(3)*sin(phi1)) ((U1(2)^2)+((U1(1)^2)+(U1(3)^2))*cos(phi1)) (U1(2)*U1(3)*(1-cos(phi1))-U1(1)*sin(phi1)); (U1(1)*U1(3)*(1-cos(phi1))-U1(2)*sin(phi1)) (U1(2)*U1(3)*(1-cos(phi1))+U1(1)*sin(phi1)) ((U1(3)^2)+((U1(1)^2)+(U1(2)^2))*cos(phi1))]; %3D Rotation Matrix

if phi1 == 0;
   XYZu = XYZt;
   Vu = Vp;
else
   XYZu = (R1*XYZt')'; %Apply Rotation Matrix to Translated Stain Locations
   Vu = (R1*Vp')'; %Apply Rotation Matrix to Projected Stain Velocity Vectors
end

Vu_test = Vu;

%Create Plane for Plotting
w = null(Sn); %Find two Orthonormal Vectors which are Orthogonal to v
[P,Q] = meshgrid(-50:1:50); %Provide a Gridwork (you choose the size)
Sx = Xp(1)+w(1,1)*P+w(1,2)*Q; %Compute the Corresponding Cartesian Coordinates using the two Vectors in w
Sy = Xp(2)+w(2,1)*P+w(2,2)*Q; %Compute the Corresponding Cartesian Coordinates using the two Vectors in w
Sz = Xp(3)+w(3,1)*P+w(3,2)*Q; %Compute the Corresponding Cartesian Coordinates using the two Vectors in w

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Create Trajectories  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating Trajectory Endpoints Into Room Dimension ONLY
if InOutTrajectory == 1
    for ij = 1:numstains
        %Trajectory Into Room Dimension ONLY
        x01(ij,:) = XYZu(ij,1);
        z01(ij,:) = XYZu(ij,3);
        x12(ij,:) = XYZu(ij,1)-1000000*Vu(ij,1); %Selecting Projected Alpha Direction
        z12(ij,:) = XYZu(ij,3)-1000000*Vu(ij,3); %Selecting Projected Alpha Direction
    end
else
    %Create Trajectory Endpoints Into & Out of Room Dimensions
    for ij = 1:numstains
        %Trajectory In/Out of Room Dimension
        x01(ij,:) = XYZu(ij,1)-1000000*Vu(ij,1); %Selecting Projected Alpha Direction
        z01(ij,:) = XYZu(ij,3)-1000000*Vu(ij,3); %Selecting Projected Alpha Direction
        x12(ij,:) = XYZu(ij,1)+1000000*Vu(ij,1); %Selecting Projected Alpha Direction
        z12(ij,:) = XYZu(ij,3)+1000000*Vu(ij,3); %Selecting Projected Alpha Direction
    end
end

numstains = sum(~isnan(alpha_p)); %Remove NaNs
alpha(isnan(x12)) = []; %Remove NaNs
minor(isnan(x12)) = []; %Remove NaNs
x01(isnan(x12)) = []; %Remove NaNs
z01(isnan(x12)) = []; %Remove NaNs
alpha_p(isnan(x12)) = [];  %Remove NaNs
alpha_pg(isnan(x12)) = []; %Remove NaNs
z12(isnan(x12)) = []; %Remove NaNs
x12(isnan(x12)) = []; %Remove NaNs

x_cs2 = x12; %Renaming variables
z_cs2 = z12; %Renaming variables
x_cs1 = x01; %Renaming variables
z_cs1 = z01; %Renaming variables

trajectory = [x_cs1 z_cs1 x_cs2 z_cs2]; %Compiling Trajectories

numstains = size(trajectory,1); %Revaluating Variables

Refp = Ref-Xp; %Translate User Defined Reference Point to Origin

if phi1 == 0;
   Refu = Refp;
else
   Refu = (R1*Refp')'; %Rotate User Defined Reference Point to XZ-plane
end

x_ref = Refu(1); %Select X-coordinate of User Defined Reference Point
y_ref = Refu(2); %Select Y-coordinate of User Defined Reference Point
z_ref = Refu(3); %Select Z-coordinate of User Defined Reference Point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Determine Bisector Intersections  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select Trajectories
xi = combvec2(x_cs1', x_cs1')'; %Compiling All Possible Trajectory Initial X-Coordinate Combinations and Purmutations
zi = combvec2(z_cs1', z_cs1')'; %Compiling All Possible Trajectory Initial Z-Coordinate Combinations and Purmutations
idb = any(diff(sort([xi zi],2),[],2)==0,2);
zi(idb,:)=[]; %Remove Repeated Combinations
xf = combvec2(x_cs2', x_cs2')'; %Compiling All Possible Trajectory Final X-Coordinate Combinations and Purmutations
xf(idb,:)=[]; %Remove Repeated Combinations
zf = combvec2(z_cs2', z_cs2')'; %Compiling All Possible Trajectory Final Z-Coordinate Combinations and Purmutations
zf(idb,:)=[]; %Remove Repeated Combinations
Alpha = combvec2(alpha',alpha')'; %Compiling All Possible Alpha Combinations and Purmutations
Minor = combvec2(minor',minor')'; %Compiling All Possible Minor Axes Combinations and Purmutations
Face = combvec2(face',face')'; %Compiling All Possible Face Combinations and Purmutations
V = combvec2(v',v')'; %Compiling All Possible Vector Combinations and Purmutations
V(idb,:)=[]; %Remove Repeated Combinations
Face(idb,:)=[]; %Remove Repeated Combinations
Alpha(idb,:)=[]; %Remove Repeated Combinations
Minor(idb,:)=[]; %Remove Repeated Combinations
xi(idb,:)=[]; %Remove Repeated Combinations
integ = [1:numstains]; %Compiling Integers
Integ = combvec2(integ,integ)'; %Compiling All Possible Trajectory Integers
Integ(idb,:)=[]; %Remove Repeated Combinations
[~,ida] = unique(sort([xi zi],2),'rows'); %Remove Repeated Purmutations
xi = xi(sort(ida),:); %Remove Repeated Purmutations
V = V(sort(ida),:); %Remove Repeated Purmutations
Face = Face(sort(ida),:); %Remove Repeated Purmutations
Alpha = Alpha(sort(ida),:); %Remove Repeated Purmutations
Minor = Minor(sort(ida),:); %Remove Repeated Purmutations
zi = zi(sort(ida),:); %Remove Repeated Purmutations
xf = xf(sort(ida),:); %Remove Repeated Purmutations
zf = zf(sort(ida),:); %Remove Repeated Purmutations
Integ = Integ(sort(ida),:); %Remove Repeated Purmutations

if nnz([xi zi]) == 0
    Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
    return
end

XZ1 = [xi(:,1) zi(:,1) xf(:,1) zf(:,1)]; %Compile Trajectories
XZ2 = [xi(:,2) zi(:,2) xf(:,2) zf(:,2)]; %Compile Adjacent Trajectories
XZ1(isinf(XZ1) | isnan(XZ1)) = []; %Remove Infinities and NaNs
XZ2(isinf(XZ2) | isnan(XZ2)) = []; %Remove Infinities and NaNs

%Iterating through all Trajectories to Determine Trajectory Intersections
for kk =1:length(xi)
    lineSegmentIntersect(XZ1(kk,:),XZ2(kk,:)); %Function Determining Line Segement Intersections: [Xc1 Zc1 Xs1 Zs1] intersecting [Xc2 Zc2 Xs2 Zs2]
    Xint12(kk) = ans.intMatrixX; %X-Coordinate Trajectory Intersection between Adjacent & Trajectory traj_shift Intersection
    Zint12(kk) = ans.intMatrixY; %Z-Coordinate Trajectory Intersection between Adjacent & Trajectory traj_shift Intersection
end

Xint12 = Xint12'; %All Possible X-Coordinate Trajectory Intersection Combinations and Permutations
Zint12 = Zint12'; %All Possible Z-Coordinate Trajectory Intersection Combinations and Permutations

xi(isnan(Xint12),:) = []; %Remove Stain Trajectories that Do Not Intersect
zi(isnan(Xint12),:) = []; %Remove Stain Trajectories that Do Not Intersect
xf(isnan(Xint12),:) = []; %Remove Stain Trajectories that Do Not Intersect
zf(isnan(Xint12),:) = []; %Remove Stain Trajectories that Do Not Intersect
Xint12(isnan(Xint12)) = []; %Remove Stain Trajectories that Do Not Intersect
Yint12 = zeros(size(Xint12));
Zint12(isnan(Zint12)) = []; %Remove Stain Trajectories that Do Not Intersect

Xc1 = xi(:,1); %Select Trajectory Endpoints
Zc1 = zi(:,1); %Select Trajectory Endpoints
Xc2 = xi(:,2); %Select Trajectory Endpoints
Zc2 = zi(:,2); %Select Trajectory Endpoints
Xs1 = xf(:,1); %Select Trajectory Endpoints
Zs1 = zf(:,1); %Select Trajectory Endpoints
Xs2 = xf(:,2); %Select Trajectory Endpoints
Zs2 = zf(:,2); %Select Trajectory Endpoints

Int12 = [Xint12, Zint12]; %Compiling Trajectory Intersection Coordinates
Int12(isnan(Int12)) = []; %Remove Stain Trajectories that Do Not Intersect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Determine Bisector Lines  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A12 = [Xc1, Zc1]; %Compiling Trajectory Endpoint Coordinates
B12 = [Xs2, Zs2]; %Compiling Trajectory Endpoint Coordinates
C12 = [Xs1, Zs1]; %Compiling Trajectory Endpoint Coordinates
D12 = [Xc2, Zc2]; %Compiling Trajectory Endpoint Coordinates

ABx = [Xc1 Xint12 Xs2 Xc1]; %Compiling Bisector Triangle Combinations for Selection
ABz = [Zc1 Zint12 Zs2 Zc1]; %Compiling Bisector Triangle Combinations for Selection
BCx = [Xs2 Xint12 Xs1 Xs2]; %Compiling Bisector Triangle Combinations for Selection
BCz = [Zs2 Zint12 Zs1 Zs2]; %Compiling Bisector Triangle Combinations for Selection
CDx = [Xs1 Xint12 Xc2 Xs1]; %Compiling Bisector Triangle Combinations for Selection
CDz = [Zs1 Zint12 Zc2 Zs1]; %Compiling Bisector Triangle Combinations for Selection
DAx = [Xc2 Xint12 Xc1 Xc2]; %Compiling Bisector Triangle Combinations for Selection
DAz = [Zc2 Zint12 Zc1 Zc2]; %Compiling Bisector Triangle Combinations for Selection

%Determine which Bisector is in Direction of Reference Point
for tri = 1:size(ABx,1)
    inAB(tri) = inpolygon(x_ref, z_ref, ABx(tri,:), ABz(tri,:)); %Determine whether Reference Point Lies within Single Bisector Triangle Combination
    inBC(tri) = inpolygon(x_ref, z_ref, BCx(tri,:), BCz(tri,:)); %Determine whether Reference Point Lies within Single Bisector Triangle Combination
    inCD(tri) = inpolygon(x_ref, z_ref, CDx(tri,:), CDz(tri,:)); %Determine whether Reference Point Lies within Single Bisector Triangle Combination
    inDA(tri) = inpolygon(x_ref, z_ref, DAx(tri,:), DAz(tri,:)); %Determine whether Reference Point Lies within Single Bisector Triangle Combination
end

%Select Bisector in Direction of Reference Point
for trj = 1:size(ABx,1)
    if inAB(trj) == 1
        endpts(trj,:) = [A12(trj,:),B12(trj,:)]; %Select Bisector from Single Bisector Triangle Combination Containing Reference Point
    elseif inBC(trj) == 1
        endpts(trj,:) = [B12(trj,:),C12(trj,:)]; %Select Bisector from Single Bisector Triangle Combination Containing Reference Point
    elseif inCD(trj) == 1
        endpts(trj,:) = [C12(trj,:),D12(trj,:)]; %Select Bisector from Single Bisector Triangle Combination Containing Reference Point
    elseif inDA(trj) == 1
        endpts(trj,:) = [D12(trj,:),A12(trj,:)]; %Select Bisector from Single Bisector Triangle Combination Containing Reference Point
    else
        Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
        ind = NaN;
        return
    end
end



if nnz(endpts(:,1)) < 2
    Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
    ind = NaN;
    return
end

numstains = size(endpts,1); %Revaluating Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Determine Trajectory Intersection Angles, \delta  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %Calculate Alpha_gravity from Intersection Points
%     t0 = [Xint12, Zint12]; %Bisector Trajectory Intersection Point
%     t1 = [endpts(:,1), endpts(:,2)]; %Endpoint of First Trajectory being Bisected
%     t2 = [endpts(:,3), endpts(:,4)]; %Endpoint of Second Trajectory being Bisected
%     nT1 = (t2 - t0) / norm(t2 - t0); %Normalize Vectors
%     nT2 = (t1 - t0) / norm(t1 - t0); %Normalize Vectors
% 
%     for ti = 1:numstains
%             delta(ti,:) = atan2(norm(det([nT2(ti,:); nT1(ti,:)])), dot(nT1(ti,:), nT2(ti,:))); %Angle between lines P0-P1 & P0-P2 
%     end

%Iterating through all Intersecting Trajectories to Determine Bisector Angle Vectors
for iii = 1:size(endpts,1)
    px12A = endpts(iii,1); %Trajectory Endpoint
    px12B = Int12(iii,1); %Trajectory Intersection Point
    px12C = endpts(iii,3); %Trajectory Endpoint
    pz12A = endpts(iii,2); %Trajectory Endpoint
    pz12B = Int12(iii,2); %Trajectory Intersection Point
    pz12C = endpts(iii,4); %Trajectory Endpoint
    v12A = [(px12A-px12B) (pz12A-pz12B)]; %Vector between Trajectory Endpoints and Trajectory Intersection Point
    v12B = [(px12C-px12B) (pz12C-pz12B)]; %Vector between Trajectory Endpoints and Trajectory Intersection Point
    v12A = v12A/norm(v12A); %Normalizing Vector
    v12B = v12B/norm(v12B); %Normalizing Vector
    bisect12pos = (v12A+v12B)*numstains^2; %Bisector Angle Vector Into the Cast-off Circle
    bisect12neg = -(v12A+v12B)*numstains^2; %Bisector Angle Vector Out of the Cast-off Circle
    px1(iii) = px12A; %Collect Trajectory Intersection Point and Endpoints
    px2(iii) = px12B; %Collect Trajectory Intersection Point and Endpoints
    px3(iii) = px12C; %Collect Trajectory Intersection Point and Endpoints
    pz1(iii) = pz12A; %Collect Trajectory Intersection Point and Endpoints
    pz2(iii) = pz12B; %Collect Trajectory Intersection Point and Endpoints
    pz3(iii) = pz12C; %Collect Trajectory Intersection Point and Endpoints
    bxpos12(iii) = bisect12pos(1)*1000000*max_room_size; %Collect Trajectory Bisectors
    bzpos12(iii) = bisect12pos(2)*1000000*max_room_size; %Collect Trajectory Bisectors
    bxneg12(iii) = bisect12neg(1)*1000000*max_room_size; %Collect Trajectory Bisectors
    bzneg12(iii) = bisect12neg(2)*1000000*max_room_size; %Collect Trajectory Bisectors
    
    bisect12A(iii,:) =[px12B pz12B] ; %Trajectory Intersection Point to Project Bisector Angle from
    bisect12B(iii,:) = bisect12pos*1000000*max_room_size; %Bisector Angle Vector Into the Cast-off Circle

%Plotting All Trajectory Intersection Arms with Corresponding Bisector
% plot([px3(iii),px2(iii)],[pz3(iii),pz2(iii)],'Color','g','LineWidth',2) %PlottingSelected Trajectories for Referencing
% plot([px1(iii),px2(iii)],[pz1(iii),pz2(iii)],'Color','g','LineWidth',2) %PlottingSelected Trajectories for Referencing
% h9(iii) = plot([px2(iii),bxpos12(iii)+px2(iii)],[pz2(iii),bzpos12(iii)+pz2(iii)],'Color','c','LineWidth',1); %Plotting Selected Trajectories Bisector Angle Vectors for Referencing (+)
% % plot([px2(iii),bxneg12(iii)+px2(iii)],[pz2(iii),bzneg12(iii)+pz2(iii)],'Color','c','LineWidth',2) %Plotting Selected Trajectories Bisector Angle Vectors for Referencing (-)
end

py2 = zeros(size(px2));
bypos12 = zeros(size(bxpos12));

%Plotting Specific, define "test" at start, Trajectory Intersection Arms with Corresponding Bisector
% plot([px3(test),px2(test)],[pz3(test),pz2(test)],'Color','g','LineWidth',2) %PlottingSelected Trajectories for Referencing
% plot([px1(test),px2(test)],[pz1(test),pz2(test)],'Color','g','LineWidth',2) %PlottingSelected Trajectories for Referencing
% plot([px2(test),bxpos12(test)+px2(test)],[pz2(test),bzpos12(test)+pz2(test)],'Color','c','LineWidth',2) %Plotting Selected Trajectories Bisector Angle Vectors for Referencing (+)
% % plot([px2(test),bxneg12(test)+px2(test)],[pz2(test),bzneg12(test)+pz2(test)],'Color','c','LineWidth',2) %Plotting Selected Trajectories Bisector Angle Vectors for Referencing (-)

bisect12 = [bisect12A (bisect12B+[px12B.*ones(numstains,1) pz12B.*ones(numstains,1)])]; %Bisector Angle Vector from and between Tradjectory & Adjacent Trajectory traj_shift Intersection
nan_row_b = any(isnan(bisect12), 2); %Remove NaN Rows (False statements) from Resuts
bisect12((nan_row_b), :) = []; %Remove NaN Rows (False statements) from Resuts

%Note: All Replaced Elements are Removed from Final Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Determine Bisector Intersection Points (Cast-off Center)  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Bisect12 = combvec2(bisect12', bisect12')'; %Compiling All Possible Bisector Intersection Combinations and Permutations
Endpts = combvec2(endpts', endpts')'; %Compiling All Possible Endpoint Combinations and Permutations
Endpts(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Px2 = combvec2(px2, px2)'; %Compiling All Possible Endpoint Combinations and Permutations
Px2(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Bxpos12 = combvec2(bxpos12, bxpos12)'; %Compiling All Possible Endpoint Combinations and Permutations
Bxpos12(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Pz2 = combvec2(pz2, pz2)'; %Compiling All Possible Endpoint Combinations and Permutations
Pz2(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Bzpos12 = combvec2(bzpos12, bzpos12)'; %Compiling All Possible Endpoint Combinations and Permutations
Bzpos12(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Xi = combvec2(xi', xi')'; %Compiling All Possible Bisector Initial X-Coordinate Combinations and Permutations
Xi(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Zi = combvec2(zi', zi')'; %Compiling All Possible Bisector Initial Z-Coordinate Combinations and Permutations
Zi(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Xf = combvec2(xf', xf')'; %Compiling All Possible Bisector Final X-Coordinate Combinations and Permutations
Xf(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Zf = combvec2(zf', zf')'; %Compiling All Possible Bisector Final Z-Coordinate Combinations and Permutations
Zf(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
INTEG = combvec2(Integ', Integ')'; %Compiling All Possible Integers
INTEG(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
Bisect12(any(diff(sort(Bisect12,2),[],2)==0,2),:)=[]; %Remove Repeated Combinations
[~,bint,~] = unique(sort(Bisect12,2),'rows'); %Remove Repeated Permutations
Bisect12 = Bisect12(sort(bint),:); %Remove Repeated Permutations
Endpts = Endpts(sort(bint),:); %Remove Repeated Permutations
Px2 = Px2(sort(bint),:); %Remove Repeated Permutations
Bxpos12 = Bxpos12(sort(bint),:); %Remove Repeated Permutations
Pz2 = Pz2(sort(bint),:); %Remove Repeated Permutations
Bzpos12 = Bzpos12(sort(bint),:); %Remove Repeated Permutations
Xi = Xi(sort(bint),:); %Remove Repeated Permutations
Zi = Zi(sort(bint),:); %Remove Repeated Permutations
Xf = Xf(sort(bint),:); %Remove Repeated Permutations
Zf = Zf(sort(bint),:); %Remove Repeated Permutations
INTEG = INTEG(sort(bint),:); %Remove Repeated Permutations

bisect1 = Bisect12(:,1:4); %Breakup Vector into Initial Components
bisect2 = Bisect12(:,5:8); %Breakup Vector into Final Components

%Iterating through all Intersecting Bisector Angles to Determine Center Coordinate Locations
for jjj = 1:size(Bisect12,1)
    lineSegmentIntersect(bisect1(jjj,:),bisect2(jjj,:)); %Function Determining the Intersection Points of Adjacent (1 to 2) & Every Other traj_shift (2 to 3) Bisecting Angles
    XintB12(jjj,:) = ans.intMatrixX; %X-Coordinate Location Results Vector from Intersection of Adjacent (1 to 2) & Every Other traj_shift (2 to 3) Bisecting Angles
    ZintB12(jjj,:) = ans.intMatrixY; %Z-Coordinate Location Results Vector from Intersection of Adjacent (1 to 2) & Every Other traj_shift (2 to 3) Bisecting Angles
end

if exist('XintB12') == 0
    Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
    ind = NaN;
    return
end

if nnz(XintB12) < 1
    Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
    ind = NaN;
    return
end

roomx = [aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)]; %Determine X-coordinates of Room Polygon
roomz = [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)]; %Determine Z-coordinates of Room Polygon

%Determine whether Bisectors Intersections within Room Polygon
if InRoom == 1 
    for in = 1:length(XintB12)
        inroom(in) = inpolygon(XintB12(in), ZintB12(in), roomx, roomz); %Select Bisector Intersections within Room Polygon
    end

    XintB12(~inroom) = []; %Remove Bisector Intersections that Lie Outside of Room Polygon
    ZintB12(~inroom) = []; %Remove Bisector Intersections that Lie Outside of Room Polygon
end

numstains = length(XintB12); %Revaluating Variable

%Cast-off Center of Origin Calculation
DD = XintB12;
c_sum_x = DD; %Compiling Center Coordinate Location Results from all Three Selected Trajectory Intersections into Single Vectornan_row = any(isnan(c_sum_x), 2); %Remove NaN Rows (False statements) from Resuts %c_sum_x = DD+rand(size(DD,1),size(DD,2)); %WITH NOISE Compiling Center Coordinate Location Results from all Three Selected Trajectory Intersections into Single Vector
nan_row_x = any(isnan(c_sum_x),2); %Remove NaN Rows (False statements) from Resuts
zero_row_x = any(c_sum_x==0,2); %Remove Zeroes (False statements) from Resuts
inf_row_x = any(c_sum_x==0987654321,2); %Remove Infinities (False statements) from Resuts
EE = ZintB12;
c_sum_z = EE; %Compiling Center Coordinate Location Results from all Three Selected Trajectory Intersections into Single Vector
nan_row_z = any(isnan(c_sum_z), 2); %Remove NaN Rows (False statements) from Resuts
zero_row_z = any(c_sum_z==0,2); %Remove Zeroes (False statements) from Resuts
inf_row_z = any(c_sum_z==0987654321,2); %Remove Infinities (False statements) from Resuts

nan_row = logical(nan_row_x.*nan_row_z);
zero_row = logical(zero_row_x.*zero_row_z);
inf_row = logical(inf_row_x.*inf_row_z);
row_check = logical(sum([nan_row,zero_row,inf_row],2));

%Remove Elements where Bisector Vectors Do Not Intersect
c_sum_x((row_check),:) = []; %Remove NaN, Zero, and Inf Rows (False statements) from Resuts
c_sum_z((row_check),:) = []; %Remove NaN, Zero, and Inf Rows (False statements) from Resuts
Xi((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
Zi((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
Xf((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
Zf((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
XintB12((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
ZintB12((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
bisect1((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
bisect2((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
INTEG((nan_row),:) = []; %Remove NaN Rows (False statements) from Resuts
INTEG = INTEG(:,1:2); %Collect First Set of Integers
BB = c_sum_x; %Renaming Variables
CC = c_sum_z; %Renaming Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Determine Cast-off Radius  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use bisector intersection point as reference point [XintB12 ZintB12]
%Use trajectory initial point as first vector point [x_cs1 z_cs1]
%Use trajectory end point as second vector point [x_cs2 z_cs2]

numstains = size(XintB12,1); %Revaluating Variable

%Iterating through all Intersecting Bisector Angles to Determine Distance from Respective Trajectories
for kkk = 1:numstains
    p2lpt1223 = [XintB12(kkk), 0, ZintB12(kkk)]; %Intersection of Bisector Angles of Every Other traj_shift Trajectories & Every Other traj_shift
    p2lva1223 = [x_cs1(kkk), 0, z_cs1(kkk)]; %Endpoint of Line between Intersections of Adjacent Trajectories & Every Other traj_shift Bisecting Angles
    p2lvb1223 = [x_cs2(kkk), 0, z_cs2(kkk)]; %Endpoint of Line between Intersections of Adjacent Trajectories & Every Other traj_shift Bisecting Angles
    point_to_line(p2lpt1223, p2lva1223, p2lvb1223); %Function Determining the Shortest Distance to a Point from a Line
    r_cv1223(kkk) = ans.r_cv; %Radius Results Vector from Intersection of Adjacent (1 to 2) & Every Other traj_shift*2 (2 to 3) Bisecting Angles 
end

r_cv1223 = r_cv1223'; %Radius Results Vector from Intersection of Adjacent (1 to 2) & (2 to 3) Bisecting Angles 
r_sum = r_cv1223; %Compiling Radius Results from all Three Selected Trajectory Intersections into Single Vector
RR = r_sum; %Renaming Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Determine Trajectory Slope & Normal Point through Center Point  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Determine Slope for each Relative Stain Trajectory
% for ijk = 1:numstains
%     slope_m(ijk) = (Zf(ijk)-Zi(ijk))/(Xf(ijk)-Xi(ijk)); %Slope of Stain Trajectories
% end
% 
% perp_m = -1.*(1./slope_m);
% 
% for sec = 1:numstains
%     point = [BB(sec) CC(sec)]; %Define Point
%     intcpt = point(2) - perp_m(sec).*point(1); %Calculate z-azis intercept
%     xvct = point(1)-Direction(1):point(1)+Direction(3); %Ã‚â€˜XÃ‚â€™ Vector for Tangent
%     tngt = perp_m(sec).*xvct + intcpt; %Calculate Z points
%     xvct_1(sec,:) = xvct;
%     tngt_1(sec,:) = tngt; %calculate Z points
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Determine Bisector Intersection Angles  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %Calculate Bisector Trajectory Intersection Angle
%     P0 = [BB, CC]; %Bisector Trajectory Intersection Point
%     P1 = [Endpt1, Endpt2]; %Endpoint of First Trajectory being Bisected
%     P2 = [Endpt3, Endpt4]; %Endpoint of Second Trajectory being Bisected
%     n1 = (P2 - P0) / norm(P2 - P0); %Normalize Vectors
%     n2 = (P1 - P0) / norm(P1 - P0); %Normalize Vectors
% 
%     for bi = 1:numstains
%             zeta(bi,:) = atan2(norm(det([n2(bi,:); n1(bi,:)])), dot(n1(bi,:), n2(bi,:))); %Angle between lines P0-P1 & P0-P2 
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Determine Distance between Bisector & Intersection  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Distance between Stain Trajectory Inntersection & Bisector Intersection Points
for bid = 1:numstains
    dist1(bid,:) = sqrt(((XintB12(bid)-bisect1(bid,1))^2)+((ZintB12(bid)-bisect1(bid,2))^2)); %Distance Equation for Distance between Stain Trajectory Inntersection & Bisector Intersection Points for Bisector #1
    dist2(bid,:) = sqrt(((XintB12(bid)-bisect2(bid,1))^2)+((ZintB12(bid)-bisect2(bid,2))^2)); %Distance Equation for Distance between Stain Trajectory Inntersection & Bisector Intersection Points for Bisector #2
end

% bsctintdist = min([dist1 dist2],[],2); %Choose the Shortest Distance between the Two Bisectors
% bsctintdist = max([dist1 dist2],[],2); %Choose the Largest Distance between the Two Bisectors
bsctintdist = mean([dist1 dist2],2); %Choose the Average Distance between the Two Bisectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Remove Potential Splashing & Large Alpha  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Choose Center Locations from Trajectories with Alpha Values between 30-60 degrees
% if alpha30less == 1 && alpha60more == 1
%     BB((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     CC((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     RR((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
% %     zeta((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     XintB12((FF),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     ZintB12((FF),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     Px2((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bxpos12((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Pz2((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bzpos12((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     BB((GG),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     CC((GG),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     RR((GG),:) = []; %Remove Locations with Alpha Values above 60 degrees
% %     zeta((GG),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     XintB12((GG),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     ZintB12((GG),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     Px2((GG),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bxpos12((GG),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Pz2((GG),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bzpos12((GG),:) = []; %Remove Bisector that Lie Outside of Room Polygon
% elseif alpha30less == 1 && alpha60more ~= 1
%     BB((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     CC((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     RR((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
% %     zeta((FF),:) = []; %Remove Locations with Alpha Values below 30 degrees
%     XintB12((FF),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     ZintB12((FF),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     Px2((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bxpos12((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Pz2((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bzpos12((FF),:) = []; %Remove Bisector that Lie Outside of Room Polygon
% elseif alpha30less ~= 1 && alpha60more == 1
%     BB((HH),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     CC((HH),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     RR((HH),:) = []; %Remove Locations with Alpha Values above 60 degrees
% %     zeta((HH),:) = []; %Remove Locations with Alpha Values above 60 degrees
%     XintB12((HH),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     ZintB12((HH),:) = []; %Remove Bisector Intesections that Lie Outside of Room Polygon
%     Px2((HH),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bxpos12((HH),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Pz2((HH),:) = []; %Remove Bisector that Lie Outside of Room Polygon
%     Bzpos12((HH),:) = []; %Remove Bisector that Lie Outside of Room Polygon
% end

d_alpha = ((((sin(alpha_orig).^2)./(1-(sin(alpha_orig).^2))).*((1+(sin(alpha_orig).^2)).*((((stdev*0.1).^2)*ones(numstains,1))./((minor*0.1).^2))))); %Alpha Impact Angle Uncertainty in Degrees
d_alpha(d_alpha>45*pi/180)=45*pi/180; %Maximum allowed Alpha Uncertainty
dalpha = sqrt(sum(d_alpha.^2)/numstains); %Root Sum Squared Alpha Impact Angle Uncertainty in Degrees
if sum(deltagamma) == 0
    d_gamma = (0.4204*exp(0.0541*alpha_orig*180/pi))*pi/180; %Gamma Impact Angle Uncertainty in Degrees
else
    d_gamma = deltagamma;
end
dgamma = sqrt(sum((d_gamma.^2))/numstains); %Root Sum Squared Gamma Impact Angle Uncertainty in Degrees
del_rad = sqrt(sum((dalpha)^2+(dgamma)^2)); %Total Uncertainty of Alpha Impact Angle

if any(d_alpha < dalpha_range(1),1) || any(d_alpha > dalpha_range(2),1)
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "d_alpha" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('dalpha =','',num2str(d_alpha),' (ideal "dalpha" Range =','',num2str(dalpha_range))); %Inputs are within Predefined Range
end
if any(d_gamma < dgamma_range(1),1) || any(d_gamma > dgamma_range(2),1)
  fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "d_gamma" is not ideal.'); %Inputs are within Predefined Range
  fprintf(fileID, '%s\r\n', strcat('dgamma =','',num2str(d_gamma),' (ideal "dgamma" Range =','',num2str(dgamma_range))); %Inputs are within Predefined Range
end

if any(del_rad) == 0;
    Psi_tot = zeros(size(isocubes)); %End Iteration if Cluster is Empty
    ind = NaN;
    clstr_grade = NaN;
    return
end

center_x = mean(BB); %Resulting Center Location X-Coordinate Vector for all Three Selected Trajectory Intersections
center_y = zeros(size(center_x)); %Selected Center Location Y-Coordinate Vector for all Three Selected Trajectory Intersections
center_z = mean(CC); %Resulting Center Location X-Coordinate Vector for all Three Selected Trajectory Intersections
center = [center_x,center_y,center_z];
radius = mean(RR); %Resulting Radius Vector for all Three Selected Trajectory Intersections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Determine Trajectory Slope Normal Point through Center Point  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine Slope for each Relative Stain Trajectory
for ijk = 1:numstains
    slope_m1(ijk,:) = (Zf(ijk,1)-Zi(ijk,1))/(Xf(ijk,1)-Xi(ijk,1)); %Slope of Stain Trajectories
    slope_m2(ijk,:) = (Zf(ijk,2)-Zi(ijk,2))/(Xf(ijk,2)-Xi(ijk,2)); %Slope of Stain Trajectories
end

b = -1;
a1 = slope_m1; %Slope of Stain Trajectories
a2 = slope_m2; %Slope of Stain Trajectories
perp_m1 = -1.*(1./slope_m1); %Perpendicular Slope of Stain Trajectories
perp_m2 = -1.*(1./slope_m2); %Perpendicular Slope of Stain Trajectories

for sec = 1:numstains
    point1 = [Xi(sec,1) Zi(sec,1)]; %Define Point
    point2 = [Xi(sec,2) Zi(sec,2)]; %Define Point
    intcpt1 = point1(2) - slope_m1(sec).*point1(1); %Calculate z-azis intercept
    intcpt2 = point2(2) - slope_m2(sec).*point2(1); %Calculate z-azis intercept
    c1(sec,:) = intcpt1; %Store intercept
    c2(sec,:) = intcpt2; %Store intercept
end

x_line1 = (b.*((b.*BB)-(a1.*CC))-(a1.*c1))./((a1.^2)+(b.^2)); %Endpoint of First Trajectory being Bisected
z_line1 = (a1.*((-b.*BB)+(a1.*CC))-(b.*c1))./((a1.^2)+(b.^2)); %Endpoint of First Trajectory being Bisected
x_line2 = (b.*((b.*BB)-(a2.*CC))-(a2.*c2))./((a2.^2)+(b.^2)); %Endpoint of Second Trajectory being Bisected
z_line2 = (a2.*((-b.*BB)+(a2.*CC))-(b.*c2))./((a2.^2)+(b.^2)); %Endpoint of Second Trajectory being Bisected
xcheck1 = x_line1-BB;
zcheck1 = z_line1-CC;
xcheck2 = x_line2-BB;
zcheck2 = z_line2-CC;

for L = 1:numstains
    if xcheck1(L) >= 0 && zcheck1(L) >= 0
        theta1(L) = atan(perp_m1(L));
    elseif xcheck1(L) < 0 && zcheck1(L) >= 0
        theta1(L) = pi+atan(perp_m1(L));
    elseif xcheck1(L) < 0 && zcheck1(L) < 0
        theta1(L) = pi+atan(perp_m1(L));
    elseif xcheck1(L) >= 0 && zcheck1(L) < 0
        theta1(L) = 2*pi+atan(perp_m1(L));
    end
    if xcheck2(L) >= 0 && zcheck2(L) >= 0
        theta2(L) = atan(perp_m2(L));
    elseif xcheck2(L) < 0 && zcheck2(L) >= 0
        theta2(L) = pi+atan(perp_m2(L));
    elseif xcheck2(L) < 0 && zcheck2(L) < 0
        theta2(L) = pi+atan(perp_m2(L));
    elseif xcheck2(L) >= 0 && zcheck2(L) < 0
        theta2(L) = 2*pi+atan(perp_m2(L));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Create and Plot Cluster Arcs  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Theta from Cast-off Center & Trajectory Points Perpendicular to Trajectoy through Center Location
s0 = [BB, CC]; %Bisector Trajectory Intersection Point
s1 = [x_line1, z_line1]; %Endpoint of First Trajectory being Bisected
s2 = [x_line2, z_line2]; %Endpoint of Second Trajectory being Bisected
ns1 = (s2 - s0) / norm(s2 - s0); %Normalize Vectors
ns2 = (s1 - s0) / norm(s1 - s0); %Normalize Vectors

for si = 1:numstains
    eta(si,:) = atan2(norm(det([ns2(si,:); ns1(si,:)])), dot(ns1(si,:), ns2(si,:))); %Angle between lines P0-P1 & P0-P2 
end

for sj = 1:numstains
    if abs(theta1(sj)-theta2(sj))-1e-6 <= eta(sj) && eta(sj) <= abs(theta1(sj)-theta2(sj))+1e-6;
        lin1(sj) = min(theta1(sj),theta2(sj)); %Initial Eta Range
        lin2(sj) = max(theta1(sj),theta2(sj)); %Final Eta Range
    else
        lin1(sj) = (min(theta1(sj),theta2(sj))-eta(sj)); %Initial Eta Range
        lin2(sj) = min(theta1(sj),theta2(sj)); %Final Eta Range
    end
end

lin12 = combvec2(lin1,lin2)';

dtheta = (0.1*res)/radius;
[M,I] = max(abs(lin12(:,1)-lin12(:,2)));

lin_1 = lin12(I,1);
lin_2 = lin12(I,2);

if (abs(lin_2-lin_1)/dtheta) < 1
    d_theta = 1;
else
    d_theta = abs(lin_2-lin_1)/dtheta;
end

if d_theta > max_room_size;
    d_theta = max_room_size;
end

for sk = 1:numstains
    theta(sk,:) = linspace(lin1(sk),lin2(sk),d_theta); %Equally Space Points between Lines
end

x_arc = center_x+radius.*cos(theta); %Arc X-coordinate
y_arc = center_y.*ones(size(x_arc)); %Arc Y-coordinate
z_arc = center_z+radius.*sin(theta); %Arc Z-coordinate

lin1 = lin_1;
lin2 = lin_2;

Eta = max([max(eta) abs(lin2-lin1)]); %Angle between Reference Lines

theta_ref = lin2+((2*pi-(lin2-lin1))/2); %Determine Reference Line furtherst Outside of Arc between lin1 and lin2
if theta_ref>(2*pi)
    theta_ref = theta_ref-(2*pi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Plot Results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
chi = linspace(0, 2*pi, 100); %Angle Vector

x_fin = (center(1) + radius.*cos(chi))'; %Resultant X-coordinate of Cast-off Circle
y_fin = (center(2)*ones(size(x_fin))); %Selected Y-coordinate of Cast-off Circle
z_fin = (center(3) + radius.*sin(chi))'; %Resultant Z-coordinate of Cast-off Circle

%Rotate Results from XZ-plane to Original Best Plane
U2 = U1; %Choose Previously Determined Line of Intersection between XZ-plane and Plane Best Fitting Stain Velocity Vectors 
phi2 = -phi1; %Planar Angle to Rotate Plane about Line of Intersection to XZ-plane
R2 = [((U2(1)^2)+((U2(2)^2)+(U2(3)^2))*cos(phi2)) (U2(1)*U2(2)*(1-cos(phi2))-U2(3)*sin(phi2)) (U2(1)*U2(3)*(1-cos(phi2))+U2(2)*sin(phi2)); (U2(1)*U2(2)*(1-cos(phi2))+U2(3)*sin(phi2)) ((U2(2)^2)+((U2(1)^2)+(U2(3)^2))*cos(phi2)) (U2(2)*U2(3)*(1-cos(phi2))-U2(1)*sin(phi2)); (U2(1)*U2(3)*(1-cos(phi2))-U2(2)*sin(phi2)) (U2(2)*U2(3)*(1-cos(phi2))+U2(1)*sin(phi2)) ((U2(3)^2)+((U2(1)^2)+(U2(2)^2))*cos(phi2))]; %3D Rotation Matrix
pu = [px2' py2' pz2'];
bupos12 = [bxpos12' bypos12' bzpos12'];
xarc = x_arc'; %X-coordinate of Arc
yarc = y_arc'; %Y-coordinate of Arc
zarc = z_arc'; %Z-coordinate of Arc
Arcu = [xarc(:) yarc(:) zarc(:)]; %Reconstructed Cast-off Arc

if phi2 == 0;
   Intt = [Xint12 Yint12 Zint12];
   Fint = [x_fin y_fin z_fin];
   centert = center';
   pt = pu;
   btpos12 = bupos12;
   Arct = Arcu;
else
   Intt = (R2*[Xint12 Yint12 Zint12]')'; %Apply Reversed Rotation Matrix to Stain Velocity Vector Intersections
   Fint = (R2*[x_fin y_fin z_fin]')'; %Apply Reversed Rotation Matrix to Cast-off Radius & Circle
   centert = R2*center'; %Apply Reversed Rotation Matrix to Cast-off Center Location
   pt = (R2*pu')';
   btpos12 = (R2*bupos12')';
   Arct = (R2*Arcu')'; %Translated Reconstructed Cast-off Arc
end

% % For Plotting Triangles to Select Triangle Containing User Defined Reference Point
% XYZnt = (R2*XYZu')';
% Vnt = (R2*Vu')';
% ABXt = (R2*(ABx(:,1:3))')';
% ABZt = (R2*(ABz(:,1:3))')';
% BCXt = (R2*(BCx(:,1:3))')';
% BCZt = (R2*(BCz(:,1:3))')';
% CDXt = (R2*(CDx(:,1:3))')';
% CDZt = (R2*(CDz(:,1:3))')';
% DAXt = (R2*(DAx(:,1:3))')';
% DAZt = (R2*(DAz(:,1:3))')';
% reft = (R2*[x_ref y_ref z_ref]')';

%Translate Stains from Origin to Original Distance from Origin
Intp = Intt+[Xp(1).*ones(size(Xint12,1),1) Xp(2).*ones(size(Xint12,1),1) Xp(3).*ones(size(Xint12,1),1)]; %Apply Reversed Translation to Stain Velocity Vector Intersections
Finp = Fint+[Xp(1).*ones(size(x_fin,1),1) Xp(2).*ones(size(x_fin,1),1) Xp(3).*ones(size(x_fin,1),1)]; %Apply Reversed Translation to Cast-off Radius & Circle
centerp = centert'+Xp; %Apply Reversed Translation to Cast-off Center Location
pp = pt+[Xp(1).*ones(size(pt,1),1) Xp(2).*ones(size(pt,1),1) Xp(3).*ones(size(pt,1),1)];
bppos12 = btpos12+[Xp(1).*ones(size(btpos12,1),1) Xp(2).*ones(size(btpos12,1),1) Xp(3).*ones(size(btpos12,1),1)];
Arcp = Arct+[Xp(1).*ones(size(Arct,1),1) Xp(2).*ones(size(Arct,1),1) Xp(3).*ones(size(Arct,1),1)];
% % For Plotting Triangles to Select Triangle Containing User Defined Reference Point
% XYZnp = XYZnt+[Xp(1).*ones(size(XYZnt,1),1) Xp(2).*ones(size(XYZnt,1),1) Xp(3).*ones(size(XYZnt,1),1)];
% ABXp = ABXt+Xp;
% ABZp = ABZt+Xp;
% BCXp = BCXt+Xp;
% BCZp = BCZt+Xp;
% CDXp = CDXt+Xp;
% CDZp = CDZt+Xp;
% DAXp = DAXt+Xp;
% DAZp = DAZt+Xp;
% refp = reft+Xp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Distribute Weight to Adjacent Cubes  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Weight = (1/dalpha).*(1/dgamma); %Weighting by Alpha, Gamma, and Zeta (prod(sum(sin(zeta).^2)))

fprintf(fileID, '%s\r\n', strcat('____')); % adds blank line for readability
fprintf(fileID, '%s\r\n', strcat('alpha degrees =',' ',num2str(round(alpha_orig(1)*180/pi)),', ', num2str(round(alpha_orig(2)*180/pi)),', ', num2str(round(alpha_orig(3)*180/pi)) )); % prints the alpha values of the 
fprintf(fileID, '%s\r\n', strcat('Weight =',' ',num2str(round(Weight))));

theta12 = [theta1' theta2'];
for intgr = 1:numstains
  INTGR = theta12(INTEG == intgr);
  INTGR = INTGR(1);
  int_stain(intgr,:) = center+radius*[cos(INTGR) 0 sin(INTGR)];
  dist_stain(intgr,:) = sqrt((int_stain(intgr,1)-XYZu(intgr,1))^2+(int_stain(intgr,2)-XYZu(intgr,2))^2+(int_stain(intgr,3)-XYZu(intgr,3))^2);
end

cosangle3=rms([dot(Sn,v(1,:)),dot(Sn,v(2,:)),dot(Sn,v(3,:))]); %Off-plane dot product
angle3=pi/2-abs(acos(cosangle3)); %Off-plane Angle
angular_delta=(rms([d_alpha,d_gamma])); %this is root mean square
Spread_Fact_cu= 2^0.5*rms([tan(angular_delta), tan(angle3)])*rms(dist_stain); %Spreading Factor %Uncertainty in Distance between a given Cube and Arc in centimeters

if Spread_Fact_cu<SF_cu_range(1)
  Spread_Fact_cu = SF_cu_range(1); %Theta Spreading Factor Floor
elseif Spread_Fact_cu>SF_cu_range(2)
  Spread_Fact_cu = SF_cu_range(2); %Theta Spreading Factor Floor
else

end
    
if Spread_Fact_cu >= SF_cu_range(1) && Spread_Fact_cu < SF_cu_range(2)
 fprintf(fileID, '%s\r\n', strcat('Spread_Fact_cu (cm)=','',num2str(round(Spread_Fact_cu)))); %Inputs are within Predefined Range
else
 fprintf(fileID, '%s\r\n', 'WARNING: The inputted value for "Spread_Fact_cu" is not ideal.'); %Inputs are within Predefined Range
 fprintf(fileID, '%s\r\n', strcat('Spread_Fact_cu (cm)=','',num2str(round(Spread_Fact_cu)),' (ideal "Spread_Fact_cu" Range =','',num2str(SF_cu_range))); %Inputs are within Predefined Range
end

for cu3 = 1:numel(cu_cx)
    cu_v(cu3,:) = [cu_cx(cu3)-centerp(1),cu_cy(cu3)-centerp(2),cu_cz(cu3)-centerp(3)]; %Vector between Cube Centers and Arc Center
    proj_cuv(cu3,:) = cu_v(cu3,:)-((cu_v(cu3,:)*Sn')./sqrt(sum(Sn.^2,2)).^2).*Sn; %Projected Vector cu_v onto Mean Plane
    
    if phi1 == 0;
       cup = [proj_cuv(cu3,:); proj_cuv(cu3,:); proj_cuv(cu3,:)];
    else
       cup = [R1*[proj_cuv(cu3,:); proj_cuv(cu3,:); proj_cuv(cu3,:)]']'; %Translate and Rotate Projected Vector to XZ-plane
    end
    
    CUP(cu3,:) = cup(1,:); %Store Translated and Rotated Projected Vector to XZ-plane
    slope_p(cu3,:) = (CUP(cu3,3)/CUP(cu3,1)); %Determine Slope of Translated and Rotated Projected Vector
    perp_p(cu3,:) = -1*(1/slope_p(cu3,:)); %Determine Perpendicular Slope of Translated and Rotated Projected Vector
    xcheckp(cu3,:) = CUP(cu3,1);
    zcheckp(cu3,:) = CUP(cu3,3);
    
    %Theta Corrections
    if xcheckp(cu3) >= 0 && zcheckp(cu3) >= 0
        thetap(cu3,:) = abs(atan(slope_p(cu3)));
    elseif xcheckp(cu3) < 0 && zcheckp(cu3) >= 0
        thetap(cu3,:) = pi-abs(atan(slope_p(cu3)));
    elseif xcheckp(cu3) < 0 && zcheckp(cu3) < 0
        thetap(cu3,:) = 0.5*3*pi-abs(atan(slope_p(cu3)));
    elseif xcheckp(cu3) >= 0 && zcheckp(cu3) < 0
        thetap(cu3,:) = 2*pi-abs(atan(slope_p(cu3)));
    end
    cu_arc(cu3,:) = centerp+radius*proj_cuv(cu3,:)/sqrt(sum(proj_cuv(cu3,:).^2,2)); %Intersection of Projected Vector and Arc (Closest Point to Arc from Cube)
    dist_cu(cu3,:) = sqrt((cu_arc(cu3,1)-cu_cx(cu3))^2+(cu_arc(cu3,2)-cu_cy(cu3))^2+(cu_arc(cu3,3)-cu_cz(cu3))^2); %Shortest Distance from Cube Centers to Arc
    CU_ARC(cu3,:) = [centerp(1)-cu_cx(cu3), centerp(2)-cu_cy(cu3), centerp(3)-cu_cz(cu3)]; %Vector between Castoff Center and Cube Center
    upsilon(cu3,:) = acos(dot(CU_ARC(cu3,:),proj_cuv(cu3,:),2)/(sqrt(sum(CU_ARC(cu3,:).^2,2))*sqrt(sum(proj_cuv(cu3,:).^2,2)))); %Off-plane Angle Upsilon
    upsilon(any(upsilon(cu3,:)>(pi),2)) = pi-upsilon(any(upsilon(cu3,:)>(0.5*pi),2));
    
    %Distance between Cube and Reconstructed Cast-off Arc PDF Distribution
    Psi_cu(cu3,:) = (1/(Spread_Fact_cu*sqrt(2*pi)))*exp((-1/2)*(dist_cu(cu3)/(Spread_Fact_cu)).^2); %Weight to be Added to Each Cube
    
end

%Product Distribution
Psi_tot = Psi_cu;

%Product Distribution Floor
Psi_tot(any(Psi_tot<(max(Psi_tot)*exp(-0.5*(sig_n^2))),2)) = max(Psi_tot)*exp(-0.5*(sig_n^2)); %PDF Floor to avoid Diminishing Spatial Regions

%Plot Room Dimensions, Stains, and Stain Trajectories
if ishandle(3)
      figure(3)
  title(['Cast-off Center & Radius of Origin: ',num2str(ip), ' of ', num2str(clstr_num(iq)),' (', num2str(iq), '/', num2str(comb_num),')'])
  view([Sn])
  h3 = quiver3(Xs,Ys,Zs,v(:,1),v(:,2),v(:,3),'Color','r');
  h4 = plot3(Xp(1),Xp(2),Xp(3),'.','MarkerSize',25,'Color','b');
  h5 = quiver3(Xp(1),Xp(2),Xp(3),Sn(1)*25,Sn(2)*25,Sn(3)*25,'Color','b');
  h6 = surf(Sx,Sy,Sz,'FaceAlpha',0.2,'EdgeAlpha',0.2); %Plot the Resulting Plane
  for ref1 = 1:numstains;
      h11(ref1) = plot3([(Xs(ref1)-10000*v(ref1,1)) (Xs(ref1)+10000*v(ref1,1))],[(Ys(ref1)-10000*v(ref1,2)) (Ys(ref1)+10000*v(ref1,2))],[(Zs(ref1)-10000*v(ref1,3)) (Zs(ref1)+10000*v(ref1,3))],'Color','r','LineWidth',1);
  end;
else
  figure(3)
  hold on
  grid on
  h1_front = plot3([aoi(1) aoi(1) aoi(1) aoi(1) aoi(1)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','c','LineWidth',5); %Plot Front Surface Dimensions
  h1_downward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(5) aoi(5) aoi(5) aoi(5) aoi(5)],'Color','g','LineWidth',4); %Plot Downward Surface Dimensions
  h1_back = plot3([aoi(2) aoi(2) aoi(2) aoi(2) aoi(2)], [aoi(4) aoi(3) aoi(3) aoi(4) aoi(4)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','b','LineWidth',3); %Plot Back Surface Dimensions
  h1_upward = plot3([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(3) aoi(3) aoi(4) aoi(4) aoi(3)], [aoi(6) aoi(6) aoi(6) aoi(6) aoi(6)],'Color','y','LineWidth',2); %Plot Downward Surface Dimensions
  h2 = plot3(Xs,Ys,Zs, '.','MarkerSize',max_room_size*0.1,'Color','r','LineWidth',2); %Plot XYZ Stains
  title(['Cast-off Center & Radius of Origin: ',num2str(ip), ' of ', num2str(clstr_num(iq)),' (', num2str(iq), '/', num2str(comb_num),')'])
  xlabel(['X-Axis (cm)'])
  ylabel(['Y-Axis (cm)'])
  zlabel(['Z-Axis (cm)'])
  view([Sn])
  h3 = quiver3(Xs,Ys,Zs,v(:,1),v(:,2),v(:,3),'Color','r');
  h4 = plot3(Xp(1),Xp(2),Xp(3),'.','MarkerSize',25,'Color','b');
  h5 = quiver3(Xp(1),Xp(2),Xp(3),Sn(1)*25,Sn(2)*25,Sn(3)*25,'Color','b');
  h6 = surf(Sx,Sy,Sz,'FaceAlpha',0.2,'EdgeAlpha',0.2); %Plot the Resulting Plane
  % h7 = plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'.','MarkerSize',25,'Color','g');
  % h8 = quiver3(XYZp(:,1),XYZp(:,2),XYZp(:,3),Vp(:,1),Vp(:,2),Vp(:,3),'Color','g');
  % h9 = plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'.','MarkerSize',25,'Color','m');
  % h10 = quiver3(XYZp(:,1),XYZp(:,2),XYZp(:,3),Vu(:,1),Vu(:,2),Vu(:,3),'Color','m');
  for ref1 = 1:numstains;
      h11(ref1) = plot3([(Xs(ref1)-10000*v(ref1,1)) (Xs(ref1)+10000*v(ref1,1))],[(Ys(ref1)-10000*v(ref1,2)) (Ys(ref1)+10000*v(ref1,2))],[(Zs(ref1)-10000*v(ref1,3)) (Zs(ref1)+10000*v(ref1,3))],'Color','r','LineWidth',1);
  end;
  % h12 = plot3(Intp(:,1),Intp(:,2),Intp(:,3),'.','MarkerSize',max_room_size*0.05,'Color','b','LineWidth',2); %Plot All Trajectory Intersection Points for Refernce
  % h13 = plot3(Finp(:,1),Finp(:,2),Finp(:,3),':r','LineWidth',3); %Plot Resultant Cast-off Circle
  % h14 = plot3(centerp(1),centerp(2),centerp(3),'*','MarkerSize',10,'Color','r','LineWidth',3); %Plot Resultant Cast-off Center Location
  % h15 = plot3(x_actual,y_actual,z_actual,'Color','m','LineWidth',3); %Plot Actual Cast-off Circle
  % h16 = plot3(actual_x,actual_y,actual_z,'p','MarkerSize',10,'Color','m','LineWidth',3); %Plot Actual Cast-off Center Location
  % for bsct = 1:size(px2,2)
  %     h17(bsct) = plot3([pp(bsct,1),bppos12(bsct,1)+pp(bsct,1)],[pp(bsct,2),bppos12(bsct,2)],[pp(bsct,3),bppos12(bsct,3)+pp(bsct,3)],'Color','k','LineWidth',1); %Plotting Selected Trajectories Bisector Angle Vectors for Referencing (+)
  % end
  % h18 = plot3(XintB12(1:end), zeros(size(XintB12)), ZintB12(1:end),'.','MarkerSize',10,'Color','g','LineWidth',3); %Plot Bisector Intersections
  % for ref1 = 1:size(XYZp)
  %     plot3([(XYZu(ref1,1)-1000*v(ref1,1)) (XYZu(ref1,1)+1000*v(ref1,1))],[(XYZu(ref1,2)-1000*v(ref1,2)) (XYZu(ref1,2)+1000*v(ref1,2))],[(XYZu(ref1,3)-1000*v(ref1,3)) (XYZu(ref1,3)+1000*v(ref1,3))],'Color','red');
  % end
  h = [h1_front h1_downward h1_back h1_upward h2 h3 h4 h5 h6]; % h7 h8 h9 h10 h11 h12 h13 h14 h15 h16 h17 h19
  legend(h, 'Front Surface', 'Downward Surface', 'Back Surface', 'Upward Surface', 'Spatter Stains', 'Stain Trajectories','Point Closest to Clustered Stains','Normal Vector Minimizing nxv','Clustered Plane', 'Location', 'northeastoutside'); % h7:'Stains Projected onto Clustered Plane', h8:'Stain Velocity Vectors Projected onto Clustered Plane', 'Stain Trajectories Projected to Clustered Plane', h9:'Clustered Stain Projections Rotated to XZ-Plane', h10:'Clustered Stain Velocity Vectors Rotated to XZ-Plane', h11:'Stain Trajectory Intersection Points', h13:'Resultant Cast-off Circle', h14:'Resultant Cast-off Center Location', h15:'Actual Cast-off Circle', h16:'Actual Cast-off Center Location', h17:'Stain Trajectory Bisector Vectors', h19:'Resultant Cast-off Arc' 
  xlim([aoi(1)-10,aoi(2)+10]);
  ylim([aoi(3)-10,aoi(4)+10]);
  zlim([aoi(5)-10,aoi(6)+10]);
  set(gca,'FontSize',20);
  axis equal;
end

Lengths = [res,res,res];
for noah = 1:numstains
  if phi1 == 0;
     Arcp = [x_arc(noah,:); y_arc(noah,:); z_arc(noah,:)]'+[Xp(1).*ones(size(x_arc,2),1) Xp(2).*ones(size(x_arc,2),1) Xp(3).*ones(size(x_arc,2),1)];
  else
     Arcp = (R2*[x_arc(noah,:); y_arc(noah,:); z_arc(noah,:)])'+[Xp(1).*ones(size(x_arc,2),1) Xp(2).*ones(size(x_arc,2),1) Xp(3).*ones(size(x_arc,2),1)];
  end
     Arc(:,:,noah) = Arcp;
        h19 = plot3(Arcp(:,1),Arcp(:,2),Arcp(:,3),'LineWidth',5,'Color','r');
    for cu = 1:size(theta,2)
        if aoi(1)<=Arcp(cu,1) && Arcp(cu,1)<=aoi(2) && aoi(3)<=Arcp(cu,2) && Arcp(cu,2)<=aoi(4) && aoi(5)<=Arcp(cu,3) && Arcp(cu,3)<=aoi(6)
            sub_col(noah,cu) = ceil((Arcp(cu,1)-xmin)*Nx/Lx);
            sub_sht(noah,cu) = ceil((Arcp(cu,2)-ymin)*Ny/Ly);
            sub_row(noah,cu) = ceil((Arcp(cu,3)-zmin)*Nz/Lz);
        else
            sub_col(noah,cu) = NaN;
            sub_sht(noah,cu) = NaN;
            sub_row(noah,cu) = NaN;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Plot Gaussian Distributions  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gaussian Distribution for Impact Angle, \alpha, Weighting
% figure(4)
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Alpha12,fa,'.')
% grid on
% title('Bell Curve')
% xlabel('Impact Angles \alpha')
% ylabel('Gauss Distribution') 

%Gaussian Distribution for Bisector Intersection Angle, \delta, Weighting
% figure(5)
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Delta,fd,'.')
% grid on
% title('Bell Curve')
% xlabel('Trajectory Intersection Angles \delta')
% ylabel('Gauss Distribution') 

%Gaussian Distribution for Stain Trajectory Intersection Angle, \zeta, Weighting
% figure(6)
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(zeta,fz,'.')
% grid on
% title('Bell Curve')
% xlabel('Bisector Intersection Angles \zeta')
% ylabel('Gauss Distribution')

% %Gaussian Distribution for Distance between Stain Trajectory Intersections & Bisector Intersections, d, Weighting
% figure(7)
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(bsctintdist,fb,'.')
% grid on
% title('Bell Curve')
% xlabel('Distance between Trajectory Intersection & Bisector Intersection')
% ylabel('Gauss Distribution')
% axis([0, 10*max_room_size, 0, max(fb)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Plot Triangle Containing User Defined Reference Point  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [min(Psi_cu) max(Psi_cu) min(dist_cu) max(dist_cu); min(Psi_theta) max(Psi_theta) min(thetap) max(thetap); min(Psi_upsilon) max(Psi_upsilon) min(upsilon) max(upsilon)];
% [del_rad*Spread_Fact_cu del_rad*Spread_Fact_theta*(lin2-lin1) del_rad*Spread_Fact_upsilon Spread_Fact_cu Spread_Fact_theta Spread_Fact_upsilon del_rad];
% 
%  figure(12); hold on; grid on
%  plot(dist_cu,Psi_cu,'.')
%  title({'\Psi_c_u_,_a_r_c v. Cube Dist for Cluster #',ip})
%  ylabel(['\Psi_c_u_,_a_r_c'])
%  xlabel(['Distance between Cube & Arc (cm)'])
%  hold off
% % 
%    figure(13); hold on; grid on
%    plot(thetap,Psi_theta,'.')
%    title({'\Psi_\theta v. \theta',ip})
%    ylabel(['\Psi_\theta'])
%    xlabel(['\theta (radians)'])
%    hold off
% % 
%  figure(14); hold on; grid on
%  plot(upsilon,Psi_upsilon,'.')
%  title({'\Psi_\Upsilon v. \Upsilon',ip})
%  ylabel(['\Psi_\Upsilon'])
%  xlabel(['Upsilon (radians)'])
%  hold off
%  
%  figure(15); hold on; grid on
%  plot(del_rad,((1/(Spread_Fact_delrad*sqrt(2*pi)))*exp((-1/2)*(del_rad/Spread_Fact_delrad)^2)),'.')
%  title({'\Psi_\delta_rad v. \delta_rad',ip})
%  ylabel(['\Psi_\delta_rad'])
%  xlabel(['\delta_rad'])
%  hold off
% 
% figure(15); hold on; grid on
% hist(Psi_cu,50)
% title({'\Psi Histogram',ip})
% xlabel(['\Psi'])
% % ylabel(['Distance between Cube & Arc (cm)'])
% hold off
% 
% figure(16); hold on; grid on
% hist(Psi_tot,100)
% title({'\Psi Histogram',ip})
% xlabel(['\Psi'])
% % ylabel(['Distance between Cube & Arc (cm)'])
% hold off

% %Plot Triangles to Select Triangle Containing User Defined Reference Point
% figure(50)
% hold on
% grid on
% axis square
% plot(ABXp(:,1),ABZp(:,3),'Color','r','LineWidth',5)
% plot(BCXp(:,1),BCZp(:,3),'Color','g','LineWidth',3)
% plot(CDXp(:,1),CDZp(:,3),'Color','b','LineWidth',2)
% plot(DAXp(:,1),DAZp(:,3),'Color','y','LineWidth',1)
% plot(refp(1,1),refp(1,3),'x','MarkerSize',20,'Color','k','LineWidth',2)
% plot(Intp(1,1),Intp(1,3),'.','MarkerSize',20,'Color','k','LineWidth',2)
% % for traj = 1:45 
% %     point1 = [x01(traj) x12(traj)]; %Create Trajectories
% %     point2 = [z01(traj) z12(traj)]; %Create Trajectories
% %     h3(traj) = plot(point1, point2); %Plot Trajectories
% % end
% plot([pp(1),bppos12(1)+pp(1)],[pp(3),bppos12(3)+pp(3)],'Color','c','LineWidth',1);
% plot([Intt(1,1) reft(1)], [Intt(1,3) reft(3)],'Color','k','LineWidth',1)
% plot([aoi(1) aoi(2) aoi(2) aoi(1) aoi(1)], [aoi(5) aoi(5) aoi(6) aoi(6) aoi(5)],'Color','k','LineWidth',2); %Plot Room XZ-Plane Dimensions
% plot(Xs,Zs, '.','MarkerSize',max_room_size*0.1,'Color','r','LineWidth',2); %Plot XZ Stains
% legend('Reference Triangle #1', 'Reference Triangle #2', 'Reference Triangle #3', 'Reference Triangle #4', 'User Defined Reference Point', 'Trajectory Intersection Point', 'Angle Bisector Vector')%,'Stain Trajectories', 'Bisector Vector')
% plot(XintB12(1), ZintB12(1),'.','MarkerSize',10,'Color','g','LineWidth',3);
% xlabel('X-Axis (cm)')
% ylabel('Z-Axis (cm)')
% axis([min([ABXp BCXp CDXp DAXp],[],'all') max([ABXp BCXp CDXp DAXp],[],'all') min([ABZp BCZp CDZp DAZp],[],'all') max([ABZp BCZp CDZp DAZp],[],'all')])
% % s = text(xmax+5,zmax,'Back Surface');
% % set(s,'Rotation',270)
% % t = text(xmin-5,zmin,'Front Surface');
% % set(t,'Rotation',90)
% % u = text(xmin+5,zmax-5,'Upward Surface');
% % v = text(xmax-85,zmin-5,'Downward Surface');
% axis square
% xlim([-100,max_room_size+100])
% ylim([-100,max_room_size+100])
% zlim([-100,max_room_size+100])
% set(gcf, 'Position', get(0, 'Screensize')); %Make Figure Full-screen
% hold off
