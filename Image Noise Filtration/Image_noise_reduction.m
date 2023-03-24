clear variables , clear clc , close all , hold on
%% Question 1
%load the image
orig = imread("car_dis.png");
doub_orig = im2double(orig);
%gray = im2gray(doub_orig);
figure(1)
fig1 = imshow(doub_orig);
%fourier
fourier = fft2(doub_orig);
logF = log(abs(fourier));
%fourier representation
figure(2)
fig2 = imshow(logF);
%shift the centre to the origin
figure(3)
fig3 = imshow(fftshift(log(abs(fourier)+1)),[]);
title('Log-amp spectrum of image "car__dis"')
xlabel('u'); ylabel('v');
print -r150 -dpng Imlogmag.png

shifted = fftshift(log(abs(fourier)+1));

%% Question 2

%uv = fourier(1,1);
%x = log(fftshift(abs(logF))+1);
%tried using a find function but it returned empty vectors...
%[row,col] = find((fftshift(abs(logF))+1)==uv);

%fftshift only mirrored the points. meaning top left (1,1) should be in the
%centre of the matrix and that is half + 1
M = 256; N = 256;
col = M/2 + 1;
row = N/2 + 1;

%let Delta be 1:

lambda = 58-49;  %the length between 2 black lines on the image
rho = round(256./lambda);


 


%% Question 3

filt = fspecial('average',[1 9]);
filtered = imfilter(doub_orig,filt);
figure(4)
fig4 = imshow(filtered);
imwrite(filtered,'IMfil.jpg')

fourier_fil = fft2(filtered);
logF_fil = log(abs(fourier_fil));
%fourier representation
figure(5)
fig5 = imshow(logF);
%shift the centre to the origin
figure(6)
fig6 = imshow(fftshift(log(abs(fourier_fil)+1)),[]);
title('Log-amp spectrum of filtered image')
xlabel('u'); ylabel('v');
print -r150 -dpng Imfillogmag.png

shifted_fil = fftshift(log(abs(fourier_fil)+1));

figure(7)
filtered2 = imfilter(doub_orig,filt,'replicate');
fig7 = imshow(filtered2);

%% Question 4

OTF = psf2otf(filt,[256 256]);
log_OTF = log(abs(OTF));
%fourier representation
figure(8)
fig8 = imshow(log_OTF);
%shift the centre to the origin
figure(9)
fig9 = imshow(fftshift(log(abs(OTF)+1)),[]);
title('Linear scale of OTF magnitude to intensity')
xlabel('u'); ylabel('v');
print -r150 -dpng OTFlogmag.png

imagi = max(abs(imag(OTF(:))));
OTF_shift = fftshift(log(abs(OTF)+1));
atten = shifted_fil((M/2)+1,(N/2)+1+rho) ./ shifted((M/2)+1,(N/2)+1+rho);

%% Question 5

% s = 9; %length of one side of rectangle around noise point
% 
% rstart = (M/2)+1-((s-1)/2);
% rend = (M/2)+1+((s-1)/2);
% cstart = (N/2)+1+rho-((s-1)/2);
% cend = (N/2)+1+rho+((s-1)/2);
% 
% H = ones(size(orig));
% H(rstart:rend,cstart:cend) = 0;

% this creates a box around only 1 point of the frequencies
%Symmetry must be maintained for other point: (N/2 + 1 - rho , M/2 + 1)

s = 27; %length of one side of rectangle around noise point


rstart1 = (M/2)+1-((s-1)/2);
rend1 = (M/2)+1+((s-1)/2);
cstart1 = (N/2)+1-rho-((s-1)/2);
cend1 = (N/2)+1-rho+((s-1)/2);


rstart2 = (M/2)+1-((s-1)/2);
rend2 = (M/2)+1+((s-1)/2);
cstart2 = (N/2)+1+rho-((s-1)/2);
cend2 = (N/2)+1+rho+((s-1)/2);
%same rstart and rend but it's fine (redundant, for my own understanding later on)
H = ones(size(orig));
H(rstart1:rend1,cstart1:cend1) = 0;
H(rstart2:rend2,cstart2:cend2) = 0;

shifted_H = fftshift(log(abs(H)+1));
applying = shifted_H .* fourier;
inverse = ifft2(applying);

imag_new = max(abs(imag(inverse(:))));
figure(10)
fig10 = imshow(inverse);
imwrite(inverse,'IMfil2.jpg')

%log-mag of new image
figure(11)
fig11 = imshow(fftshift(log(abs(fft2(inverse))+1)),[]);
title('Log-amp spectrum of image after convolution of filter H')
xlabel('u'); ylabel('v');
print -r150 -dpng Imfil2logmag.png

%log-mag of new OTF
OTF2 = psf2otf(shifted_H,[256 256]);
figure(12)
fig12 = imshow(fftshift(log(abs(OTF2)+1)),[]);
title('magnitude of OTF after convolution filter H')
xlabel('u'); ylabel('v');
print -r150 -dpng OTF2logmag.png


