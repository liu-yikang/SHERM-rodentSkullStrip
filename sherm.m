function sherm(input_filename, output_filename, descriptor_filename, animal, isotropic)

%SHERM Main function of SHERM (SHape descriptor selected Extremal Regions
%after Morphologically filtering (SHERM)
%   Inputs: 
%   input_filename: filename of the image to process
%   output_filename: name of the output file
%   descriptor_filename: filename of shape descriptors
%   animal: 'rat' or 'mouse'
%   isotropic (1 or 0): if the input image has isotropic spatial resolution

nii = load_untouch_nii(input_filename);
drs = min(nii.hdr.dime.pixdim(2:4));

% radii of open/close kernels
open_rad = 0.2:drs:0.7;
close_rad = 0.2:drs:0.5;

% set range of brain volume (mm3) and shape descriptor
desc_mat = load(descriptor_filename);
if strcmp(animal, 'rat')
    brain_vol_range = [1500, 1900];
    desc_tmp = desc_mat.polar_template_rat;
elseif strcmp(animal, 'mouse')
    brain_vol_range = [300, 550];
    desc_tmp = desc_mat.polar_template_mouse;
else
    disp('Animal name not defined.');
end

% return MSERs and shape descriptor distances
if isotropic == 0
    [desc_dist,masks] = ...
        sherm_parallel_2d(nii,descriptor_filename,drs,brain_vol_range,...
        open_rad,close_rad,animal);
else
    [desc_dist,masks] = ...
        sherm_parallel_3d(nii,descriptor_filename,drs,brain_vol_range,...
        open_rad,close_rad,animal);
end

desc_dist_min = min(desc_dist(desc_dist>0))
[ids_chan, ids_mser] = find(desc_dist >0 & desc_dist < (desc_dist_min[0] + 0.05*sum(abs(desc_tmp))));
mask = zeros(size(masks{1}(:,:,:,1)));
for i_candidate = 1:length(ids_chan)
    mask = mask + masks{ids_chan(i_candidate)}(:,:,:,ids_mser(i_candidate));
end
mask = mask > 0.1*length(ids_chan);

nii.img = double(mask);
save_untouch_nii(nii,output_filename);

end
