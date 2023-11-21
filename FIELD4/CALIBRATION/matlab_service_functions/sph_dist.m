%
% Calcola la distanza sferica tra due punti.   by Lorenzo Chiesi
%d = sph_dist(theta1,phi1,theta2,phi2)

function d = sph_dist(theta1,phi1,theta2,phi2)

d = acos( cos((pi/2)-phi1).*cos((pi/2)-phi2) + sin((pi/2)-phi1).*sin((pi/2)-phi2).*cos(theta2-theta1) );
end

