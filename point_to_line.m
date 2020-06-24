function dist = point_to_line(pt, v1, v2)
%Function Determining the Shortest Distance to a Point from a Line
      a = v1 - v2; %Subtract Vectors
      b = pt - v2; %Subtract Vectors
      d = norm(cross(a,b)) / norm(a); %Normalize Vectors and take the Cross Product
      dist.r_cv = d; %Result Extraction
end