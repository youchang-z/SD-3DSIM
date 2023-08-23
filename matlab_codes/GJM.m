function [refracted_vector] = GJM(phi, theta, vector)
rot = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
bend = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
rot_back = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
refracted_vector = rot_back*bend*rot*vector;
end

%% GJM simulates how the objective lens refract light rays  