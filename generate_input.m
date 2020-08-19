clearvars;
clf;
%path and name of csv
OutputFileName = 'dummy.xlsx';
filename = fullfile(OutputFileName);
if exist(filename, 'file')==2
  delete(filename);
end

%user inputs

size_stain=10;%average width of stain (mm)
rc=[100 150 150]'; %%define a circle,  coordinate of center circle,  cm
dc=80;%diameter of circle, cm 
ttheta=27*pi/180; %angle theta from degree to radian
nc=[0 cos(ttheta) sin(ttheta)]';%normal vector to circle (change this to change the plane of rotation of the sweeping motion of the castoff
% the classical vector is nc=[0 1 0]'
number_point_along_circle=20; % is the number of desired stains 
room_size=[250 500 250]'% define room size
% corner of room (min x, y, z)
room_corner=[0,0,0]'
%end user input

segment_length=5*max(room_size);
z_vector=[0 0 1]';
director_vector=dc/2*cross(z_vector,nc)

%plot circle
start_circle_angle=10*rand
thetas=[start_circle_angle:floor((360-start_circle_angle)/number_point_along_circle):360]'*pi/180; %angular coordinates on circle
theta2=[1:5:360];
for i=1:length(theta2)
    R = rotate_3D(director_vector, 'any', theta2(i), nc);
    circle=rc+R;
    figure(25)
    plot3(circle(1),circle(2),circle(3),'ro'); hold on
    view(-nc)
end
view(-nc)


for i=1:length(thetas)
    R = rotate_3D(director_vector, 'any', thetas(i), nc);
    long_tangent=(cross(R,nc));
    tangent=long_tangent/norm(long_tangent);
    segment_origin=rc+R;
    segment_end=segment_origin+tangent*segment_length;
    figure(25)
    % plot3(segment_origin(1),segment_origin(2),segment_origin(3),'go'); hold on
    % plot3(segment_end(1),segment_end(2),segment_end(3),'g+'); hold on
    plot3([segment_origin(1) segment_end(1)]',[segment_origin(2) segment_end(2)]',[segment_origin(3) segment_end(3)]','-o'); hold on
    hold on
    view(nc)
    
    %determination of intersections
    distance=Inf;
    alpha=NaN;
    gamma=NaN;
    counter=0;
    % intersection with top plane
    n=[0 0 1]; t=[-1 0 0];proj_tan=[tangent(1),tangent(2),0];proj_tan=proj_tan/norm(proj_tan);
    point=[0 0 room_corner(3)+room_size(3)];
    [I,check]=plane_line_intersect(n,point,segment_origin',segment_end');
    
    if check==1
        
        new_distance=norm(segment_origin'-I);
        new_alpha=90-acosd(dot(n,tangent));
        %new_gamma=-acosd(dot(-t,proj_tan));
        new_gamma=atan2(dot((cross(-t,proj_tan)),n), dot(-t,proj_tan))*180/pi;
        
         
        if new_distance<distance
            stain_x(i)=I(1);stain_y(i)=I(2);stain_z(i)=I(3);
            counter=counter+1;
            distance=new_distance;
            alpha_vec(i)=new_alpha;
            gamma_vec(i)=new_gamma;
        end
    end
    % intersection with right plane
    n=[1 0 0];    t=[0 0 1]; proj_tan=[0,tangent(2),tangent(3)];proj_tan=proj_tan/norm(proj_tan);
    point=[room_corner(1)+room_size(1) 0 0];
    [I,check]=plane_line_intersect(n,point,segment_origin',segment_end');
    if check==1
        
        new_distance=norm(segment_origin'-I);
        new_alpha=90-acosd(dot(n,tangent));
        %new_gamma=-acosd(dot(-t,proj_tan));
        new_gamma=atan2(dot((cross(-t,proj_tan)),n), dot(-t,proj_tan))*180/pi;
        if new_distance<distance
            stain_x(i)=I(1);stain_y(i)=I(2);stain_z(i)=I(3);
            counter=counter+1;
            distance=new_distance;
            alpha_vec(i)=new_alpha;
            gamma_vec(i)=new_gamma;
        end
    end
    % intersection with left plane
    n=[-1 0 0];t=[0 0 1];proj_tan=[0,tangent(2),tangent(3)];proj_tan=proj_tan/norm(proj_tan);
    point=[0 0 room_corner(1)];
    [I,check]=plane_line_intersect(n,point,segment_origin',segment_end');
    if check==1
        
        new_distance=norm(segment_origin'-I);
        new_alpha=90-acosd(dot(n,tangent));
        %new_gamma=-acosd(dot(-t,proj_tan));
        new_gamma=atan2(dot((cross(-t,proj_tan)),n), dot(-t,proj_tan))*180/pi;
        
        if new_distance<distance
            stain_x(i)=I(1);stain_y(i)=I(2);stain_z(i)=I(3);
            counter=counter+1;
            distance=new_distance;
            alpha_vec(i)=new_alpha;
            gamma_vec(i)=new_gamma;
        end
    end
    % intersection with bottom plane
    n=[0 0 -1];t=[-1 0 0]; proj_tan=[tangent(1),tangent(2),0];proj_tan=proj_tan/norm(proj_tan);
    point=[0 0 room_corner(3)];
    [I,check]=plane_line_intersect(n,point,segment_origin',segment_end');
    if check==1
        
        new_distance=norm(segment_origin'-I);
        new_alpha=90-acosd(dot(n,tangent));
        %new_gamma=-acosd(dot(-t,proj_tan));
        new_gamma=atan2(dot((cross(-t,proj_tan)),n), dot(-t,proj_tan))*180/pi;
        
        if new_distance<distance
            stain_x(i)=I(1);stain_y(i)=I(2);stain_z(i)=I(3);
            counter=counter+1;
            distance=new_distance;
            alpha_vec(i)=new_alpha;
            gamma_vec(i)=new_gamma;
        end
    end
    
    if counter>1.1
        warning('more than one intersection found')
    end
    if counter<0.9
        warning('no intersection found')
    end
    
    %stain_size
    
        width_stain(i)=size_stain;
        length_stain(i) = width_stain(i)/sind(alpha_vec(i)); %in mm
    
end
xlabel('x');ylabel('y');zlabel('z');


%save csv file
number_cells=length(theta2);
row_begin=5%where numbers begin
row_end=5+number_cells-1;%where number end
writecell({'room length'},filename,'Range','A1');
writecell({'room width'},filename,'Range','A2');
writecell({'room height'},filename,'Range','A3');
cellReference = sprintf('B%d:B%d', 1, 3); % Create cell reference
writematrix(room_size,filename,'Range',cellReference);
writecell({'Project',	'Pattern'	'Stain'	'Location x'	'Location y'	'Location z'	'Surface'	'M axis (mm)'	'm axis (mm)'	'alpha'	'beta'	'gamma'	'Intersect x'	'Intersect y'	'Intersect z'},filename,'Range','A4:O4');
cellReference = sprintf('D%d:D%d', row_begin, row_end); % Create cell reference
writematrix(stain_x',filename,'Range',cellReference);
cellReference = sprintf('E%d:E%d', row_begin, row_end); % Create cell reference
writematrix(stain_y',filename,'Range',cellReference);
cellReference = sprintf('F%d:F%d', row_begin, row_end); % Create cell reference
writematrix(stain_z',filename,'Range',cellReference);
cellReference = sprintf('j%d:j%d', row_begin, row_end); % Create cell reference
writematrix(alpha_vec',filename,'Range',cellReference);
cellReference = sprintf('L%d:L%d', row_begin, row_end); % Create cell reference
writematrix(gamma_vec',filename,'Range',cellReference);
cellReference = sprintf('H%d:H%d', row_begin, row_end); % Create cell reference
writematrix(length_stain',filename,'Range',cellReference);
cellReference = sprintf('I%d:I%d', row_begin, row_end); % Create cell reference
writematrix(width_stain',filename,'Range',cellReference);

figure(25);view(-nc)




