%function [MOVINGREG] = XYR_registerImages(MOVING,FIXED,out,pre,iter)
clear;
iter=30;out='/data/disk2/xiongyirong/test/';pre='standard_CC';
MOVING='/data/disk2/xiongyirong/template_file/MNICC07mm_XXmask_final_final.nii';
iterations=50;
% locate template CC slice 
%F = load_untouch_nii(FIXED);
%FId=find(sum(sum(F.img,2),3));     FIMG = double(squeeze(F.img(FId,:,:))); 
FIMG=double(zeros(311,260));
% for x=1:312
%     for y=1:260
%         if (x-155)*(x-155)+(y-130)*(y-130)<10000
%             FIMG(x,y)=1;
%         end
%     end
% end
width=60;height=30;
wStart=155-width/2;wEnd=154+width/2;
hStart=130-height/2;hEnd=129+height/2;
FIMG(wStart:wEnd,hStart:hEnd)=double(ones(width,height));
% locate native CC slice 
M = load_untouch_nii(MOVING);
MId=find(sum(sum(M.img,2),3));    MIMG = double(squeeze(M.img(MId,:,:)));
fixedRefObj = imref2d(size(FIMG));
movingRefObj = imref2d(size(MIMG));
% Intensity-based registration
[optimizer, metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 1.00000e-06;
optimizer.MinimumStepLength = 1.0000e-4;
optimizer.MaximumStepLength = 1.000e-2;
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.500000;

% Align centers
[xFixed,yFixed] = meshgrid(1:size(FIMG,2),1:size(FIMG,1));
[xMoving,yMoving] = meshgrid(1:size(MIMG,2),1:size(MIMG,1));
sumFixedIntensity = sum(FIMG(:));
sumMovingIntensity = sum(MIMG(:));
fixedXCOM = (fixedRefObj.PixelExtentInWorldX .* (sum(xFixed(:).*double(FIMG(:))) ./ sumFixedIntensity)) + fixedRefObj.XWorldLimits(1);
fixedYCOM = (fixedRefObj.PixelExtentInWorldY .* (sum(yFixed(:).*double(FIMG(:))) ./ sumFixedIntensity)) + fixedRefObj.YWorldLimits(1);
movingXCOM = (movingRefObj.PixelExtentInWorldX .* (sum(xMoving(:).*double(MIMG(:))) ./ sumMovingIntensity)) + movingRefObj.XWorldLimits(1);
movingYCOM = (movingRefObj.PixelExtentInWorldY .* (sum(yMoving(:).*double(MIMG(:))) ./ sumMovingIntensity)) + movingRefObj.YWorldLimits(1);
translationX = fixedXCOM - movingXCOM; %��ͼ������Ķ��� 
translationY = fixedYCOM - movingYCOM;

% Coarse alignment
initTform = affine2d();
initTform.T(3,1:2) = [translationX, translationY];

% Apply Gaussian blur
fixedInit = imgaussfilt(FIMG,1.000000);
movingInit = imgaussfilt(MIMG,1.000000);

% Normalize images
movingInit = mat2gray(movingInit);
fixedInit = mat2gray(fixedInit);    

tform = imregtform(movingInit,movingRefObj,fixedInit,fixedRefObj,'affine',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
MOVINGREG.Transformation = tform;
[MOVINGREG.RegisteredImage] = imwarp(MIMG, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
%[RegisteredImageO] = imwarp(MIMG, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
RegisteredImageO = MOVINGREG.RegisteredImage;

% Nonlinear registration
%  [MOVINGREG.DisplacementField,MOVINGREG.RegisteredImage] = imregdemons(MOVINGREG.RegisteredImage,FIMG,200,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);
 [D{1},REIMG] = imregdemons(MOVINGREG.RegisteredImage,FIMG,iterations,'AccumulatedFieldSmoothing',1,'PyramidLevels',3,'DisplayWaitbar',0);
 MOVINGREG.RegisteredImage = REIMG;
 OUTIMG{1} = REIMG;
for iiter = 2:iter
    [D{iiter},OUTIMG{iiter}] = imregdemons(MOVINGREG.RegisteredImage,FIMG,iterations,'AccumulatedFieldSmoothing',1,'PyramidLevels',3,'DisplayWaitbar',0);
    out1 = images.geotrans.internal.applyDisplacementField(MOVINGREG.RegisteredImage,D{iiter},'linear',0,0);
    MOVINGREG.RegisteredImage = out1;
end
DD = D{1};
for j = 2:size(D,2)
    %DD = DD+D{j};
    %DD = imwarp(DD,D{j});
    DD = images.geotrans.internal.applyDisplacementField(DD,D{j},'linear',0,0);
    %DD = images.geotrans.internal.applyDisplacementField(DD,D{j},'linear',0,0);
end
MOVINGREG.DisplacementField = DD;
outs1 = imwarp(RegisteredImageO,DD);

%outs1 = imwarp(RegisteredImageO,DD,'OutputView', fixedRefObj, 'SmoothEdges', true);

RIMG = outs1;MOVINGREG.fixedRefObj = fixedRefObj;
MOVINGREG.movingRefObj = movingRefObj;
%t=regexp(pre,'_','split');t=t{2};
%for i =1:4
%     subplot(2,4,i)
%     RIMG = MOVINGREG.RegisteredImage;RIMG(RIMG<0.3+i/10)=0;RIMG(RIMG>=0.3+i/10)=1;
%     Diff = sum(sum((RIMG-FIMG).^2));
%     imshowpair(FIMG,RIMG,'falsecolor','ColorChannels','green-magenta');
%     title({['ID: ',t];['thr: ',num2str(0.3+i/10),'; Diff: ',num2str(Diff)]})
%    subplot(1,4,i)
%    RIMG1 =outs1;RIMG1(RIMG1<0.3+i/10)=0;RIMG1(RIMG1>=0.3+i/10)=1;
%    Diff = sum(sum((RIMG1-FIMG).^2));
%    imshowpair(FIMG,RIMG1,'falsecolor','ColorChannels','green-magenta');
%    title({['ID: ',t];['thr: ',num2str(0.3+i/10),'; Diff: ',num2str(Diff)]})   
%end
MOVINGREG.JDNonlinear = ZCX_jacobian(MOVINGREG.DisplacementField);

% [MIMG1] = imwarp(MIMG, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
% B = imwarp(MIMG1,MOVINGREG.DisplacementField);B(B<0.5)=0; B(B>0.5)=1;
% figure
% imshowpair(FIMG(105:240,95:150),B(105:240,95:150),'falsecolor','ColorChannels','green-magenta');
%saveas(gcf,[out,filesep,pre,'_QC.png']);
%MOVINGREG.FId = FId;MOVINGREG.MId = MId;
%MOVINGREG.Fsform = [ F.hdr.hist.srow_x; F.hdr.hist.srow_y; F.hdr.hist.srow_z];
%MOVINGREG.Msform = [ M.hdr.hist.srow_x;M.hdr.hist.srow_y; M.hdr.hist.srow_z];

save([out,filesep,pre,'_Deformation'],'MOVINGREG');
%O = M;
%O.img(1,:,:) = MOVINGREG.RegisteredImage;

%O.fileprefix = [out,filesep,pre,'_Reg2d'];
%save_untouch_nii(O,O.fileprefix)
%end



%% 
function [vx,vy] = compose(ax,ay,bx,by)

    [x,y] = ndgrid(0:(size(ax,1)-1), 0:(size(ax,2)-1)); % coordinate image
    x_prime = x + ax; % updated x values
    y_prime = y + ay; % updated y values
    
    % Interpolate vector field b at position brought by vector field a
    bxp = interpn(x,y,bx,x_prime,y_prime,'linear',0); % interpolated bx values at x+a(x)
    byp = interpn(x,y,by,x_prime,y_prime,'linear',0); % interpolated bx values at x+a(x)

    % Compose
    vx = ax + bxp;
    vy = ay + byp;
    
end
%% Jacobian
function det_J = ZCX_jacobian(Displacement)

    % Gradients
    [gx_y,gx_x] = gradient(Displacement(:,:,1));
    [gy_y,gy_x] = gradient(Displacement(:,:,2));
    
    % Add identity
    gx_x = gx_x + 1;  % zero displacement should yield a transformation T = Identity (points keep their positions)
    gy_y = gy_y + 1;  % adding identity matrix here
    
    % Determinant
    det_J = gx_x.*gy_y - ...
            gy_x.*gx_y;
end
%end