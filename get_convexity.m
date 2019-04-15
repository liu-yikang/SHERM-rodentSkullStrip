function convexity = get_convexity(mask, dx, dy, dz)

ind = find(mask);
dim = size(mask);
[y,x,z] = ind2sub(dim,ind);
x = x*dy;
y = y*dx;
z = z*dz;
[~,V] = convhull(x,y,z);
convexity = sum(mask(:)*dx*dy*dz)/V;

end

