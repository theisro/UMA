%Fully reimplemented with basical trigonometric function ready to port to C++
function [tr,pr] = sph_rotation2(t,p,ts,ps)

r = 1;

%Convert to rectangular coords
x = r*cos(t)*cos(p);
y = r*sin(t)*cos(p);
z = r*sin(p);

%Rotate
xr =   x*cos(ps)*cos(ts) + y*cos(ps)*sin(ts) + z*sin(ps);
yr = - x*sin(ts)         + y*cos(ts);
zr = - x*sin(ps)*cos(ts) - y*sin(ps)*sin(ts) + z*cos(ps);

%Return to spherical coords
rr = sqrt(xr^2+yr^2+zr^2);
if ((xr == 0)&&(yr==0))
    tr = 0;
else    
    if (xr > 0)
        tr = atan(yr/xr);
    else
        tr = pi - atan(yr/xr);
    end
end
pr = asin(zr/rr);

end
