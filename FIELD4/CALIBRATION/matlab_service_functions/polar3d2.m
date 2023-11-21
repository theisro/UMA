function polar3d2(R,AxesLimit)
%
%R is a mtrix (SurfacePoint x [theta,phi,rho(complex)])
%POLAR3D(R)
%POLAR3D(R,Max)

if exist('AxesLimit','var')
    rmax = AxesLimit;
else
    rmax = squeeze(max(abs(squeeze(R(:,3)))));
end


if (size(R,1) < 500)
    %Create Theta/Phi mesh with 5° resolution
    [T,P] = meshgrid( linspace(0,2*pi,360/5) , linspace(-pi/2,pi/2,180/5));
else
    %Create Theta/Phi mesh with 1° resolution
    [T,P] = meshgrid( linspace(0,2*pi,360/1) , linspace(-pi/2,pi/2,180/1));
end

%Complete matrix wrapping around on theta dimension
R2 = R;
R2(:,1) = R2(:,1)+(2*pi);
R = [R ; R2];

%Interpolate sparse point finding magnitude and phase value on grid
Magnitude = griddata(R(:,1),R(:,2),abs(R(:,3)),T,P,'cubic');
Phase = griddata(R(:,1),R(:,2),abs(angle(R(:,3))),T,P);

%Convert to cartesian
[X,Y,Z] = sph2cart(T,P,Magnitude); 

surf(X,Y,Z,Phase,'FaceColor','interp','FaceLighting','phong');
colormap(flipdim(jet,1));
caxis([0,pi]);

xlim([-rmax, rmax]);
ylim([-rmax, rmax]);
zlim([-rmax, rmax]);

axis vis3d

end
