%function [SR_image] = Main_SR
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function runs a complete super-resolution algorithm based on the
% continuous moments analysis. The set of images is 'LR_Tiger_xx.tif' where
% xx is a number between 01 and 40. Each image is 64 x 64 pixels. The
% output super-resolved image 'SR_image' is 512 x 512 pixels.
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************

% import files
figure
for fileIndex = 1:40
    % import file
    fileNameStr = sprintf('LR_Tiger_%02i.tif', fileIndex);
    importfile(fileNameStr)
    % plot file
    subplot(5,8,fileIndex)
    imshow(fileNameStr)
end

 % import file
importfile('HR_Tiger_01.tif')
% plot file
figure
imshow(HR_Tiger_01)

%% Compute red layer of LR images

lrImage1_red = zeros(size(LR_Tiger_01), 'uint8');
lrImage1_red(:,:,1) = LR_Tiger_01(:,:,1);

figure
imshow(lrImage1_red)













% Register images
[Tx_RGB Ty_RGB]= ImageRegistration;

% Compute super-resolved image
[SR_image]= ImageFusion(Tx_RGB, Ty_RGB);