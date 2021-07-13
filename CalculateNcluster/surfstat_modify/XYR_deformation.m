clear;clc;
MIMG=im2double(imread('/data/disk2/xiongyirong/cc1.png'));
FIMG=im2double(imread('/data/disk2/xiongyirong/fix1.png')); 
%MOVING='/data/disk2/xiongyirong/template_file/MNICC07mm_XXmask_final_final.nii';
%M = load_untouch_nii(MOVING);
%MId=find(sum(sum(M.img,2),3));    MIMG = double(squeeze(M.img(MId,:,:)));

hFig      = vis.mfigure;
hFig.Name = 'MIMG';
subplot(1,5,1), imshow(MIMG,[]); title(['imageOri' ]);%fixedRefObj = imref2d(size(FIMG));
%movingRefObj = imref2d(size(MIMG));
for i=1:4
    load(['/data/disk2/xiongyirong/testdeformation/TX' num2str(i) '.mat']);
    D(:,:,2)=Tx;
    load(['/data/disk2/xiongyirong/testdeformation/TY' num2str(i) '.mat']);
    D(:,:,1)=Ty;
    DD{i}=D;
    MIMG=imwarp(MIMG,D);
    subplot(1,5,i+1), imshow(MIMG,[]); title(['image ' num2str(i)]);
    G{5-i}=xyr_inverse_DD(D);
    fprintf('========%d Deformation Inversed===========\n',i);
end

hFig      = vis.mfigure;
hFig.Name = 'FIMG inversion1';
FIMG=im2double(imread('/data/disk2/xiongyirong/fix1.png')); 
subplot(1,5,1), imshow(FIMG,[]); title(['imageOri' ]);%fixedRefObj = imref2d(size(FIMG));

for j=1:4
    for x=1:4
        for y=1:4
            for z=1:4
               
hFig      = vis.mfigure;
hFig.Name = 'FIMG inversion';
FIMG=im2double(imread('/data/disk2/xiongyirong/fix1.png')); 
subplot(1,5,1), imshow(FIMG,[]); title(['imageOri' ])               
    FIMG=imwarp(FIMG,G{1}{x});
    subplot(1,5,2), imshow(FIMG,[]); title(['image ' num2str(x)]);
    FIMG=imwarp(FIMG,G{2}{y});
    subplot(1,5,3), imshow(FIMG,[]); title(['image ' num2str(y)]);
    FIMG=imwarp(FIMG,G{3}{z});
    subplot(1,5,4), imshow(FIMG,[]); title(['image ' num2str(z)]);
    FIMG=imwarp(FIMG,G{4}{j});
    subplot(1,5,5), imshow(FIMG,[]); title(['image ' num2str(j)]);
    fname=['/data/disk2/xiongyirong/testdeformation/PNG/' num2str(x) '_' num2str(y) '_' num2str(z) '_' num2str(j) '.png'];
    saveas(gcf, fname);
    close all;
               
            end
        end
    end
end


