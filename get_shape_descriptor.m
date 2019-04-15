function [desc] = get_shape_descriptor(mask, dx, dy, dz)

% mask = bwareaopen(mask,1000);
% n_vxl = sum(mask(:));
ind = find(mask);
dim = size(mask);
[y,x,z] = ind2sub(dim,ind);
x = x*dy;
y = y*dx;
z = z*dz;
c = [mean(x), mean(y), mean(z)];
coeff = pca([x - c(1), y - c(2), z - c(3)]);

dim_idx = [];
for ii = 1:3
    dim_idx(ii) = find(coeff(:,ii) == max(coeff(:,ii)));
end

mat = eye(4);
mat(1:3, dim_idx) = coeff;

tform = affine3d(mat);

ref = imref3d(size(mask),[-c(1), dim(1)*dx-c(1)],[-c(2), dim(2)*dy-c(2)],...
            [-c(3), dim(3)*dz-c(3)]);
ref2 = imref3d(size(mask).*[1,dy/dx,dz/dx],[-c(1), dim(1)*dx-c(1)],[-c(2), dim(2)*dy-c(2)],...
            [-c(3), dim(3)*dz-c(3)]);
mask = imwarp(double(mask), ref, tform, 'OutputView', ref2)>0.5;

area_xy = sum(sum((sum(mask,3)>0)));
% mask = mask & ~imerode(mask, return3dStrel(2,1,1,1));
ind = find(mask);
dim = size(mask);
[x,y,z] = ind2sub(dim,ind);
x = x*dx;
y = y*dx;
z = z*dx;
c = [mean(x), mean(y), mean(z)];
x = x-c(1);
y = y-c(2);

% x = x-c(1);
% y = y-c(2);
% z = z-c(3);
% coord = mat*[y,x,z,ones(size(x))]';
% y = coord(1,:) - mean(coord(1,:));
% x = coord(2,:) - mean(coord(2,:));
% % z = coord(3,:) - mean(coord(3,:));
% % x = x(z>min(z)+1 & z<max(z)-1);
% % y = y(z>min(z)+1 & z<max(z)-1);

[theta,rho] = cart2pol(x,y);
N_theta = 20;
N_rho = 10;
desc = zeros(N_theta, N_rho);
for i = 1:N_theta
    for j = 1:N_rho
        desc(i,j) = sum(theta > -pi+(i-1)*2*pi/N_theta & ...
            theta < -pi+i*2*pi/N_theta & ...
            rho > (j-1)*max(rho(:))/N_rho & ...
            rho < j*max(rho(:))/N_rho);
    end
end
%desc = desc(:)*dx*dy*dz/(sum(sum((sum(mask,3)>0)))*dx*dy)^(3/2);
desc = desc(:)/area_xy^(3/2);


% mask_close = imclose(mask, return3dStrel(sum(mask(:)*dx*dy*dz).^(1/3)/10, dx, dy, dz));
% convexity = sum(mask(:))/sum(mask_close(:));
% [~,V] = convhull(x,y,z);
% convexity = sum(mask(:)*dx*dy*dz)/V;

end

