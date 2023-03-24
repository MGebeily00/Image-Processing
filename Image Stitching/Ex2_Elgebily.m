clear variables, clear clc, close all
%% Question 1
images = imageSet('imfolder');
im1 = read(images,1);
im2 = read(images,2);
im3 = read(images,3);
im4 = read(images,4);
%% Question 3
% for n = 1:4
%     for m = n+1:4
%         [p{n,m},p{m,n}] = cpselect(read(images,n),read(images,m),'Wait',true);
%     end
% end
% save psets p

load psets.mat
%% Question 4
tform2 = estimateGeometricTransform(p{2,1},p{1,2},'projective');
warp2 = imwarp(im2,tform2);
figure(1)
imshow(warp2);
imwrite(warp2,'warp2.jpg')
%% Question 5
% sz = size(im1);
tform1 = projective2d;
% %gray = im2gray(warp2); %for me to see matrix (size was too big for matlab in RGB)
% [xLimitsOut1,yLimitsOut1] = outputLimits(tform1,[1 sz(2)],[1 sz(1)]);
% [xLimitsOut2,yLimitsOut2] = outputLimits(tform2,[1 sz(2)],[1 sz(1)]);
% % X Limits
% [xlims{1,2}] = xLimitsOut1;
% [xlims{2,1}] = xLimitsOut2;
% % Y limits
% [ylims{1,2}] = yLimitsOut1;
% [ylims{2,1}] = yLimitsOut2;

RA = imref2d(size(im1));
warp1 = imwarp(im1,RA,tform1);
[new_warp2,RB] = imwarp(im2,RA,tform2);
XWorldLimits = [min(min(RA.XWorldLimits,RB.XWorldLimits)),max(max(RA.XWorldLimits,RB.XWorldLimits))];
YWorldLimits = [min(min(RA.YWorldLimits,RB.YWorldLimits)),max(max(RA.YWorldLimits,RB.YWorldLimits))];
imagesize = [fix(YWorldLimits(2)-YWorldLimits(1)) , fix(XWorldLimits(2)-XWorldLimits(1))];
%creating imref

imref = imref2d(imagesize,XWorldLimits,YWorldLimits);

%stitch
A = imwarp(im1,tform1,'OutputView',imref); 
B = imwarp(im2,tform2,'OutputView',imref); 

%alphablend = vision.AlphaBlender;
%figure(2);
imstitch12 = max(A,B);
stitch12 = imtool(imstitch12);
% stitch12 = alphablend(A,B);
% imshow(stitch12);

%% Question 6
% 4 tform objects tranforming from 1to1 (identity same as 1to0) , 2to1, 3to2, 4to3 (n to n-1)
tform(1) = projective2d;
tform(2) = estimateGeometricTransform(p{2,1},p{1,2},'projective');
tform(3) = estimateGeometricTransform(p{3,2},p{2,3},'projective');
tform(4) = estimateGeometricTransform(p{4,3},p{3,4},'projective');

tform11 = tform(1);
tform21 = tform(2);
tform31 = projective2d(tform(3).T * tform(2).T);
tform41 = projective2d(tform(4).T * tform31(1).T);
%creating tform1 (1Tn array)
tform1 = [tform11 tform21 tform31 tform41];

%% Question 7

% sz = size(im1);
% [xLimitsOut1,yLimitsOut1] = outputLimits(tform1(1),[1 sz(2)],[1 sz(1)]);
% [xLimitsOut2,yLimitsOut2] = outputLimits(tform1(2),[1 sz(2)],[1 sz(1)]);
% [xLimitsOut3,yLimitsOut3] = outputLimits(tform1(3),[1 sz(2)],[1 sz(1)]);
% [xLimitsOut4,yLimitsOut4] = outputLimits(tform1(4),[1 sz(2)],[1 sz(1)]);


RA = imref2d(size(im1));
[newwarp2,RB] = imwarp(im2,RA,tform1(2));
[newwarp3,RC] = imwarp(im3,RA,tform1(3));
[newwarp4,RD] = imwarp(im4,RA,tform1(4));
R_xlims = [RA.XWorldLimits RB.XWorldLimits RC.XWorldLimits RD.XWorldLimits];
R_ylims = [RA.YWorldLimits RB.YWorldLimits RC.YWorldLimits RD.YWorldLimits];

XWorldLimits_new = [min(R_xlims),max(R_xlims)];
YWorldLimits_new = [min(R_ylims),max(R_ylims)];
newimagesize = [fix(YWorldLimits_new(2)-YWorldLimits_new(1)) , fix(XWorldLimits_new(2)-XWorldLimits_new(1))];
%creating imref

imref_new = imref2d(newimagesize,XWorldLimits_new,YWorldLimits_new);

%stitch (code too cumbersome. try something else. use for loop...)
% A = imwarp(im1,tform1(1),'OutputView',imref_new); 
% B = imwarp(im2,tform1(2),'OutputView',imref_new); 
% C = imwarp(im3,tform1(3),'OutputView',imref_new); 
% D = imwarp(im4,tform1(4),'OutputView',imref_new); 

% Make and empty image to save the stitch
stitch1234 = zeros([newimagesize 3], 'like', im1);

% blender
alphablend = vision.AlphaBlender('Operation', 'Binary mask','MaskSource', 'Input port');

for i = 1:4  % result for non for loop was better
    
    Im = readimage(images, i);   
   
    % Transform the images
    warpedImage = imwarp(Im, tform1(i), 'OutputView', imref_new);
     
    % stitching
    stitching = imwarp(true(size(Im,1),size(Im,2)), tform1(i), 'OutputView',imref_new);
    
    % Place the images
    stitch1234 = step(alphablend, stitch1234, warpedImage, stitching);
end
% imstitch12 = max(A,B);
% imstitch34 = max(C,D);
% imstitch1234 = max(imstitch12,imstitch34);
% stitch1234 = imshow(imstitch1234);
%stitch1234 = imtool(imstitch1234);
figure(2)
imshow(stitch1234);
%imwrite(stitch1234,'4stitched.jpg')

%% Question 8

% moving points 4 to plane 1
pts4=p{4,3};
[x4to3,y4to3] = transformPointsForward(tform(4),pts4(:,1),pts4(:,2));
[x4to2,y4to2] = transformPointsForward(tform(3),x4to3,y4to3);
[x4to1,y4to1] = transformPointsForward(tform(2),x4to2,y4to2);

tform4to1=estimateGeometricTransform(p{4,3},[x4to1,y4to1],'projective');

% transofrm 3 to 1
pts3=p{3,2};

[x3to2,y3to2] = transformPointsForward(tform(3),pts3(:,1),pts3(:,2));
[x3to1,y3to1] = transformPointsForward(tform(2),x3to2,y3to2);

tform3to1=estimateGeometricTransform(p{3,2},[x3to1,y3to1],'projective');

tformto1 = [tform(1) tform(2) tform3to1 tform4to1];
% getting error E_nm by for-looping all points to get e_nm then sum it and
% square*1/k under square root (4 images means k = 4 so 1/k = 0.25)
for n=1:4
    for m=n+1:4
        e=0;
        for i=1:4
            e=e+(norm(transformPointsForward(tformto1(n),p{n,m})-...
                transformPointsForward(tformto1(m),p{m,n})))^2;
            E(n,m)=sqrt(0.25*e);
        end
    end
end

% Getting RMS:
%initialise RMS = 0 (counter)
RMS=0;
for n=1:4
    for m=n+1:4
        RMS=sqrt((1/6)*sum(sum(E.^2)));
    end
end 

%% Question 9

% the latitude given by the image description
latitude=54;             % Indicated as 54deg N on the image
y_dist=15*60;            % 15 degree in y direction
x_dist=20*60*cosd(54);   % 20 degree in x direction

%positions in pixel number, tl, tr, bl, br
actual_corners = [74 131; 1704 134; 233 2184; 1558 2193];
mapped_corners = [0 0; round(x_dist) 0; 0 round(y_dist);...
    round(x_dist) round(y_dist)];

tformGeo = fitgeotrans(actual_corners,mapped_corners,'projective');  %create transform
immapped = imwarp(stitch1234, tformGeo);    %transform the image

%crop the image
crop_rect = [131 74 700 880];
imcropped = imcrop(immapped, crop_rect);   
imshow(imcropped); 
%imwrite(imcropped, 'croppedmap.jpg');

%% Question 12

alpha1 = tformto1(2).T(1:8);
alpha2 = tformto1(3).T(1:8);
alpha3 = tformto1(4).T(1:8);

Beta = [alpha1' alpha2' alpha3']';

Alpha1 = reshape([Beta(1,:) 1],3,3);
Alpha2 = reshape([Beta(1,:) 1],3,3);
Alpha3 = reshape([Beta(1,:) 1],3,3);

tforms(1) = tform(1);
tforms(2) = projective2d(Alpha1);
tforms(3) = projective2d(Alpha2);
tforms(4) = projective2d(Alpha3);


[RMS2,E2] = bundle_adjustment(tforms,p);



