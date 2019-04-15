function [ se, b ] = return3dStrel( r,dx,dy,dz )

width = ceil(r./[dx,dy,dz]);
[x,y,z] = meshgrid(-width(1):width(1),-width(2):width(2),-width(3):width(3));
m = sqrt((x*dx).^2 + (y*dy).^2 + (z.*dz).^2);
b = m <= r;
se = strel('arbitrary',b);

end

