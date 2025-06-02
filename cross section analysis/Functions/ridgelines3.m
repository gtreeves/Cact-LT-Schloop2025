function [bgm2,fgm] = ridgelines3(U,filename,outidx)

%%% Getting the Watershed Lines (Ridge Lines) to demarcate the different
%%% areas

U1 = imtophat(U,strel('disk',20));%figure , imshow(uint8(I));
%imshow(uint8(I))
U2 = gaussFiltDU(U1);
BW = mat2gray(U);
%figure, imshow(uint8(BW));

%figure , imshow((mat2gray(U2,[70 80])).*mask1(:,1:2595));
%figure , imshow(imerode((mat2gray(U2,[70 80])).*mask1(:,1:2595)),se);

N = 10;
BW = [BW(:,end-(2*N-1):end),BW,BW(:,1:2*N)];

% BW = BW.*mask;
% figure, imshow(uint8(BW));

BW = imopen(BW,strel('disk',20));
% BW = imopen(BW,strel('disk',5)); % not good
% BW = imopen(BW,strel('disk',25)); % almost same
fgm = imregionalmax(BW);


%D = bwdist(fgm);
D = bwdist(fgm);
DL = watershed(D);
bgm2 = DL == 0;
%figure, imshow(bgm2), title('Watershed Ridge Lines)');

% U2ws = watershed(1-mDU(U2o));
% bgm2 = U2ws == 0;

end
