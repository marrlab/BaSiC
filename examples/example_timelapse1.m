clear; clc; close all;
addpath(genpath('dctool'));

%% read in time-lapse movie
images_dir = 'BaSiCPluginDemo/Demoexamples/Timelapse_brightfield/Uncorrected/';
files =dir([images_dir '*.png']);
for i = 1:length(files)  
    IF(:,:,i) = imread([images_dir files(i).name]); % original image
end
   
%% estimate flatfield and darkfield
% For brightfield images, darkfield estimation is not necessary (default)
% also, we do not need to set regularisation parameters, lambda and
% lambda_dark, but use default setted ones
flatfield = BaSiC(IF);  
basefluor =  BaSiC_basefluor(IF,flatfield);

%% plot estimated shading files
figure; subplot(121); imagesc(flatfield);colorbar; title('Estimated flatfield');
subplot(122); plot(basefluor);xlabel('Timepoint');title('Estimated mean background signal');

%% image correction
IF_corr = zeros(size(IF));
for i = 1:length(files)
    IF_corr(:,:,i) = double(IF(:,:,i))./flatfield - basefluor(i);
end

%% view original and corrected image pair
figure; subplot(121); imagesc(IF(:,:,1));colormap('gray');title('Original image');
subplot(122);imagesc(IF_corr(:,:,1));colormap('gray');title('BaSiC corrected image');
