%[theta_out, phi_out] = sph_rotation(theta_in, phi_in, theta_shift, phi_shift)
%
%Sperical rotation "TILT THEN PAN" for direction defined in spherical coords
%Work on single direction and on vector of directions

function [theta_out, phi_out] = sph_rotation(theta_in, phi_in, theta_shift, phi_shift)


%Force input coords as column vector
in_size = size(theta_in);
theta_in = theta_in(:);
phi_in = phi_in(:);

%Compute rotation matrix for rectangular coords
R1= [  cos(phi_shift)     0     sin(phi_shift)
       0                  1     0
      -sin(phi_shift)     0     cos(phi_shift) ];
  
R2= [ cos(theta_shift)  sin(theta_shift)    0
      -sin(theta_shift) cos(theta_shift)    0
      0                 0                   1 ];
R=R1*R2;

%Convert to rectangular coords
[x,y,z] = sph2cart(theta_in,phi_in,1);

%Apply rotation
xyz =  R * [x,y,z]';

%Compute back spherical coords
[theta_out,phi_out,rho] = cart2sph(xyz(1,:),xyz(2,:),xyz(3,:));
theta_out = reshape(theta_out,in_size);
phi_out = reshape(phi_out,in_size);

end
