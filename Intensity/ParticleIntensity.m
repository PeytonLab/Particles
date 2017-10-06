%ParticleIntensity 
%Author: David Podorefsky, October6th2017

%Z range to measure
zmin = 15;
zmax = 52;

[filename, path] = uigetfile('.tif','*.*');
cd(path) 
info = imfinfo(filename);
%zmax = numel(info); only if zmin = 1

B = cell(zmax-zmin+1,1); Int = cell(zmax-zmin+1,1);
L = cell(zmax-zmin+1,1); V = cell(zmax-zmin+1,1); 

for z = zmin:zmax % Mark spots to measure intensities on image
    I = imread(filename,z); %Read image
    %figure; imagesc(I)
	%I = rgb2gray(I); %Truecolor to grayscale
	Int{z-zmin+1} = I; %Layer intensity
	%G = fspecial('gaussian',[5 5],2); I = imfilter(I,G); %Apply Gaussian filter
    BW = im2bw(I,0.6); %Binarize
    %figure; imagesc(BW)
    [b,l]= bwboundaries(BW,'noholes'); %Trace particles
    L{z-zmin+1} = double(l); %Particle locations
    %figure; imagesc(l)
    B{z-zmin+1} = zeros(1,size(b,1)); %Particle intensity storage space
end

for s = 1:2
    for z = 1:zmax-zmin+1
        for n = 1:length(B{z})		
			M = Int{z}(L{z}==n); %intensities of particle pixels
			M(M==0) = []; %remove zeros			
            if s == 1 || z == 1
                B{z}(n) = sum(M); %Sum of pixel intensity of particle
            else             
               N = L{z-1}(L{z}==n); %hits of pixels in same location in previous layer
			   N(N==0) = []; %remove misses			   
               if length(N) > 0.75*length(M)
				  B{z}(n) = NaN; %mark intensity blank for current layer             
				  nn = mode(N); %particle number in previous layer
                    if  ~isnan(B{z-1}(nn)); %not connected to a third layer
                         B{z-1}(nn) = B{z-1}(nn)+sum(M-p(1)*(z-1)); %add intensity to particle in layer below
                      else %connected by more than 2 layers 
                        for zz = 1:zmax-zmin+1
						  N = L{z-zz-1}(L{z}==n);
						  N(N==0) = [];
                            if isempty(N) %cannot be traced
                                break
                            else    
                                nnn = mode(N);				      
                                if ~isnan(B{z-zz-1}(nnn))
                                    break
                                end
                            end
                        end                     
                      if ~isempty(N)
                         B{z-zz-1}(nnn) = B{z-zz-1}(nnn)+sum(M-p(1)*(z-1));
                      end
                    end                       
               else %not attached to layer below
                    B{z}(n) = sum(M-p(1)*(z-1));
               end
  
            end
        end
    end    
    if s == 1
        for z = 1:length(B)
            zint = B{z}(~isnan(B{z}));
            V{z} = mean(zint);
        end
    	p = polyfit(1:zmax-zmin+1,cell2mat(V)',1);     
    end
end

%Output results
k = 1;
for i = 1:length(B)
    ii = length(B{i});
    for j = 1:ii
        P(k) = B{i}(j);
        k = k+1;
    end
end
P(isnan(P)) = []; %Remove 
P(P>510); %Noise removal
MeanParticleIntensity = mean(P)
TotalParticleIntensity = sum(P)