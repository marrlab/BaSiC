clear; clc; close all;
addpath(genpath('dctool'));

%% read in time-lapse movie
images_dir = 'BaSiCPluginDemo/Demoexamples/Timelapse_nanog/Uncorrected/';
files =dir([images_dir '*.tif']);
for i = 1:length(files)  
    IF(:,:,i) = imread([images_dir files(i).name]); % original image
end
   
%% estimate flatfield and darkfield
[flatfield, darkfield] = BaSiC(IF,'darkfield','true','lambda',2.0,'lambda_dark',2.0);
basefluor =  BaSiC_basefluor(IF,flatfield,darkfield);

%% plot estimated shading files
figure; subplot(131); imagesc(flatfield);colorbar; title('Estimated flatfield');
subplot(132); imagesc(darkfield);colorbar; title('Estimated darkfield');
subplot(133); plot(basefluor);xlabel('Timepoint');title('Estimated background fluorescence');
%% image correction
IF_corr = zeros(size(IF));
for i = 1:length(files)
    IF_corr(:,:,i) = (double(IF(:,:,i))-darkfield)./flatfield - basefluor(i);
end

%% view original and corrected image pair
figure; subplot(121); imagesc(IF(:,:,1));colormap('gray');title('Original image');
subplot(122);imagesc(IF_corr(:,:,1));colormap('gray');title('BaSiC corrected image');