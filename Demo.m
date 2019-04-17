root_dir = strrep(mfilename('fullpath'),'example','');

descriptor = fullfile(root_dir, 'rodent_brain_polar_descriptor.mat');
input_dir = fullfile(root_dir, 'data', 'img_N3_intensity_normalized');
output_dir = fullfile(root_dir, 'data', 'brain_mask_autoext');
img_list = dir(fullfile(input_dir, '*.nii'));

animal = {'mouse', 'mouse', 'mouse', 'mouse', 'rat', 'rat', 'rat', 'rat'};
isotropic = [1,1,0,0,0,0,0,0];

for i = 1:length(img_list)
    sherm(fullfile(input_dir, img_list(i).name), fullfile(output_dir, img_list(i).name), ...
        descriptor, animal{i}, isotropic(i));
end
