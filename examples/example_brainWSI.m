clear; clc; close all;
addpath(genpath('dctool'));

%% read in time-lapse movie
images_dir = 'BaSiCPluginDemo/Demoexamples/WSI_Brain/Uncorrected_tiles/';
files =dir([images_dir '*.tif']);
for i = 1:length(files)  
    IF(:,:,i) = imread([images_dir files(i).name]); % original image
end
   
%% estimate flatfield and darkfield
% For fluorescence images, darkfield estimation is often necessary (set
% 'darkfield' to be true)
[flatfield,darkfield] = BaSiC(IF,'darkfield','true');  

%% plot estimated shading files
% note here darkfield does not only include the contribution of microscope
% offset, but also scattering light that entering to the light path, which
% adds to every image tile
figure; subplot(121); imagesc(flatfield);colorbar; title('Estimated flatfield');
subplot(122); imagesc(darkfield);colorbar;title('Estimated darkfield');

%% image correction
IF_corr = zeros(size(IF));
for i = 1:length(files)
    IF_corr(:,:,i) = (double(IF(:,:,i))-darkfield)./flatfield;
end

%% save images for stitching
save_dir = '/Users/tying84/Documents/Humboldt_research/paper/NatureSubjournal/BaSiCsoftware/BaSiCPluginDemo/Demoexamples/WSI_Brain/MatCorrected_tiles/';
for i = 1:length(files)
    imwrite(uint16(IF_corr(:,:,i)),[save_dir files(i).name]);
end
