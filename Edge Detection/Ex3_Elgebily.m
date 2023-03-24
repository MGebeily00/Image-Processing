%% name: Mostafa Elgebily
clear variables, clear clc, close all
%% Question 1
sigma = 5;                  % set the width of the Gaussian
L = 2*ceil(sigma*4)+1;     % fill in a constant to define the matrix size (4)
xymax = (L-1)/2;            % the maximum of the x and the y coordinate
xrange = -xymax:xymax;      % the range of x values
yrange = xrange;            % the range of y values
figure(1)
h = fspecial('gaussian', L, sigma);      % create the PSF matrix
C = cat(3, ones(size(h)),ones(size(h)),zeros(size(h)));
        % create a RGB matrix to define the colour of the surface plot
hd =surf(xrange,yrange,h,C,'FaceColor','interp','Facelight','phong');
                        % create the surface plot of the gaussian
camlight right          % add a light at the right side of the scene
xlim([-0.5*L 0.5*L]);   % set appropriate axis limits
ylim([-0.5*L 0.5*L]);
xlabel('x');            % add axis labels
ylabel('y');
zlabel('h(x,y)');
print -r150 -dpng ex3_1.png % print the result to file

im = imread('ang2.png');          % read the image from file
sigma = [1 5 10 20];            % fill array with sigma values
figure(2)
for i=1:4                       % do 4 times:
h = fspecial('gaussian', L, sigma(i));              % create the PSF
imfiltered = imfilter(im,h,'same','replicate');    % apply the filter
subplot(2,2,i);                 % define a subplot
imshow(imfiltered,[]);          % show the image
title(['\sigma = ' num2str(sigma(i))]);       % include a plot title
end
print -r150 -dpng filtered_1_5_10_20.png     % print the result to file


%% Question 2
sigma = 3;
%h = exp(-(x.^2 + y.^2)/(2*sigma^2))/(2*pi*sigma^2);
%hy = (y.*exp(-(x.^2 + y.^2)/(2*sigma^2)))/(2*pi*sigma^4);

L = 2*ceil(sigma*4)+1;      % fill in a constant to define the matrix size 
N = (L-1)/2;                % get the size of half of the full range
[x,y]=meshgrid(-N:N,-N:N);  % create the coordinates of a 2D orthogonal grid
xrange = -floor(length(x)/2):floor(length(x)/2);           % the range of x values
yrange = xrange;           % the range of y values
figure(3)
h = exp(-(x.^2 + y.^2)/(2*sigma^2))/(2*pi*sigma^2);         % new PSF
hy = (y.*exp(-(x.^2 + y.^2)/(2*sigma^2)))/(2*pi*sigma^4);   % derivative y
C = cat(3, ones(size(hy)),ones(size(hy)),zeros(size(hy)));
        % create a RGB matrix to define the colour of the surface plot
hyd =surf(xrange,yrange,hy,C,'FaceColor','interp','Facelight','phong');
                        % create the surface plot of the gaussian
camlight right          % add a light at the right side of the scene
xlim([-0.5*L 0.5*L]);   % set appropriate axis limits
ylim([-0.5*L 0.5*L]);
xlabel('x');            % add axis labels
ylabel('y');
zlabel('h(x,y)');
print -r150 -dpng ex3_2b.png % print the result to file

HY = psf2otf(hy,size(im));
HY_R = max(real(HY(:)));        % almost 0
HY_I = max(imag(HY(:)));        % around 0.2


HYI = imag(HY);
shiftedHYI = fftshift(HYI);
N_OTF = (size(im)-1)/2;
[X,Y]=meshgrid(-N_OTF:N_OTF,-N_OTF:N_OTF);   % create the coordinates of a 2D orthogonal grid
figure(4)
CHY = cat(3, ones(size(shiftedHYI)),ones(size(shiftedHYI)),zeros(size(shiftedHYI)));
HD =surf(X,Y,shiftedHYI,CHY,'FaceColor','interp','Facelight','phong','edgecolor','none');
                        % create the surface plot of the gaussian
camlight right          % add a light at the right side of the scene
xlim([-0.5*150 0.5*150]); % set appropriate axis limits
ylim([-0.5*150 0.5*150]);
xlabel('u');            % add axis labels
ylabel('v');
zlabel('h(u,v)');
print -r150 -dpng ex3_2d.png % print the result to file


% 2g)
imfft = fft2(im);
conv = ifft2(imfft.*HY);
matim = mat2gray(conv);
figure(5);
imshow(matim,[])
imwrite(matim,'ex3_2g.jpg')

%% Question 3
sigma = 4;

% fx
figure(6)
fx = ut_gauss(im,sigma,1,0);
fx_gray = mat2gray(fx);
imshow(fx_gray, []);
title('fx');
imwrite(fx_gray, 'ex3_3_fx.jpg');

% fy
figure(7)
fy = ut_gauss(im,sigma,0,1);
fy_gray = mat2gray(fy);
imshow(fy_gray, []);
title('fy');
imwrite(fy_gray, 'ex3_3_fy.jpg');

% fxx
figure(8)
fxx = ut_gauss(im,sigma,2,0);
fxx_gray = mat2gray(fxx);
imshow(fxx_gray, []);
title('fxx');
imwrite(fxx_gray, 'ex3_3_fxx.jpg');

% fyy
figure(9)
fyy = ut_gauss(im,sigma,0,2);
fyy_gray = mat2gray(fyy);
imshow(fyy_gray, []);
title('fyy');
imwrite(fyy_gray, 'ex3_3_fyy.jpg');

% fxy
figure(10)
fxy = ut_gauss(im,sigma,1,1);
fxy_gray = mat2gray(fxy);
imshow(fxy_gray, []);
title('fxy');
imwrite(fxy_gray, 'ex3_3_fxy.jpg');

%% Question 4

% gradient
figure(11);
gradmag = abs(sqrt(fx.^2 + fy.^2));
gradmag_gray = mat2gray(gradmag);
imshow(gradmag_gray, []);
title('f gradient magnitude');
imwrite(gradmag_gray, 'ex3_4_f_gradmag.jpg');

% laplacian
figure(12);
laplacian = fxx + fyy;
laplacian_gray = mat2gray(laplacian);
imshow(laplacian_gray, []);
title('f laplacian');
imwrite(laplacian_gray, 'ex3_4_f_laplacian.jpg');

%% Question 5

% figure(13)
% if laplacian(x>0)
%     x = 1;
% else
%     x = 0;
% end
% imshow(laplacian);

% create a binary out of laplacian with anything above 0 (positive) = 1
figure(13);
im_binary = imbinarize(laplacian, 0);      
imshow(im_binary);
imwrite(im_binary, 'ex3_5a_binary.jpg');

% 5b

figure(14)
% Structuring element (4 neighbour diamond = diameter 4 = radius 2)
SE = strel('diamond',2);
eroded = imerode(im_binary,SE);
subtracted = im_binary - eroded;
imshow(subtracted)
imwrite(subtracted,'ex3_5b_0-crossing.jpg')

%% Question 6

% 6a
figure(15)
masked = gradmag.*subtracted;
masked_gray = mat2gray(masked);
marr_mask = imshow(masked_gray, []);

% 6b
figure(16)
% imcontrast(marr_mask)
T = (0.2853 + 0.2814)/2;
masked_thresh = masked_gray > T;
imshow(masked_thresh,[]);
imwrite(masked_thresh,'ex3_6b_masked-threshold.jpg')

% 6c
figure(17)
marker = masked_gray > 0.55;
imshow(marker)

% 6d
figure(18)
mask = masked_gray > 0.1;
imshow(mask)

% 6e
figure(19)
reconstruct = imreconstruct(marker,mask);
imshow(reconstruct)
imwrite(reconstruct,'ex3_6e_reconstruct.jpg')


%% Question 7

sigma = 3;

% Canny edge:

figure(20);
canny = ut_edge(im,'s',sigma,'c','h',[0.1 0.0033]);
imshow(canny,[]);
imwrite(canny,'ex3_7_canny.jpg')

% Marr-Hildreth edge:

figure(21);
marr = ut_edge(im, 's', sigma, 'm', 'h', [0.1 0.0033]);
imshow(marr,[]);
imwrite(marr,'ex3_7_marr.jpg')
