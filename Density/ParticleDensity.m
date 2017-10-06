%Name: ParticleDensity.m, Author: David Podorefsky
%Counts particles, computes density, measures areas and circularity
%Requires importing >1 RGB tif images
[filename,path] = uigetfile('.tif','*.*','MultiSelect','on');
cd(path)
px = 3.1; %px/um
p = cell(length(filename),5); 
p{1,1} = 'Filename'; p{1,2} = 'Count'; p{1,3} = 'Density (/mm^2)'; 
p{1,4} = 'Size (um^2)'; p{1,5} = 'Circularity';
output = inputdlg('Enter save name:','Save results',[1 50]);
for i = 1:length(filename)
    p{i+1,1} = filename{i};
    I = imread(filename{i});
    I = rgb2gray(I);
    G = fspecial('gaussian',[5 5],2); I = imfilter(I,G);
    BW = im2bw(I,graythresh(I));
    [B,L] = bwboundaries(BW,'noholes');
    p{i+1,2} = size(B,1);
    p{i+1,3} = p{i+1,2}/(size(I,1)*size(I,2)/(px*1000)^2);
    stats = regionprops(BW,'area','perimeter');
    p{i+1,4} = mean([stats.Area])/px^2; 
    circularity = (([stats.Perimeter]+pi).^2)./(4*pi*[stats.Area]); 
    p{i+1,5} = mean(circularity);
end
xlswrite(output{1},p)