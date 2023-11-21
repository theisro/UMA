function polar3d(R,AxesLimit)
%
%POLAR3D(R)                     by Lorenzo Chiesi
%POLAR3D(R,Max)
%
%Example: polar3d(ones(7,12));
%

phi_div = size(R,1)-1;
theta_div = size(R,2);

theta_span = 2*pi / theta_div;
phi_span = pi / phi_div;

theta_row = (0:theta_span:2*pi);
theta = ones(phi_div+1,1)*theta_row;

phi_col = (-pi/2:phi_span:pi/2)';
phi = phi_col*ones(1,theta_div+1);

   
plot_radius = abs([R , R(:,1)]);
    
if exist('AxesLimit','var')
    rmax = AxesLimit;
else
    rmax = squeeze(max(max(plot_radius)));
end
    


phase = double(abs(angle([R , R(:,1)])));

[X,Y,Z] = sph2cart(theta,phi,plot_radius); 

surf(X,Y,Z,phase,'FaceColor','interp','FaceLighting','phong');
colormap(flipdim(jet,1));
caxis([0,pi]);

xlim([-rmax, rmax]);
ylim([-rmax, rmax]);
zlim([-rmax, rmax]);

axis vis3d

end
