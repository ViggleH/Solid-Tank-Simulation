close all
clear all 
SimSettings
% load may16_frickeZoom7a.mat
% load may25_FlexyWaterZ6-2.mat
% load oct23zoom-redoJune12.mat
% load oct23zoom-2redux.mat %ideal with RI (manufactured, wrong stabilizer)
% load june14_FlexyRIfluid5.mat %ideal with RI
% load june21_FlexyRIfluid14.mat
% load june29_frickeAtFlexy1.mat
load june30_fricke7.mat %fricke

dia = diaset;
% dia = 104; % 104mm - diameter of bore hole
wallT = 0.125*25.4; %0.125in wall - container wall thickness
gapT = (103.6-dia)/2; % gap thickness i.e. bath thickness
% d1 = 183.5+(dia/2); % distance from front edge of acrylic block to center of bore
% d2 = 344-d1; % distance from center of bore hole to detectors - 344mm = total length of acrylic block
doff = 0; % vertical offset of laser source
r0 = 1; %defines the initial ray intesity to be one

angles = linspace(span1,span2,r); % angle of inclination of each ray from laser source
r1 = (dia/2)*scale - wallT*scale; %gel radius
r2 = (dia/2)*scale;  %container radius
r3 = (dia/2)*scale + gapT*scale; %bore radius


% path = '/home/festesio/Documents/MATLAB/Merged';
% file = '/test.mat';
% load (strcat(path,file))
quality = zeros(size(mag,1),size(mag,2),size(mag,3),size(mag,4),size(mag,5),6);


parfor i = 1:numSims
    profiles(i,:) = profile{i};
    scaleFactor = 12000/mean(profiles(i,:));
    profiles(i,:) = profiles(i,:)*scaleFactor;
    
    L2(i) = norm(profiles(i,:) - mean(profiles(i,:)),2)/(320)^0.5;
end


% quality(:,:,1) = 1./dyRange;
% quality(:,:,:,:,:,1) = -0.5*tanh(0.2*(noZbeamUni10 - 5)*2*pi) + 0.5;
% quality(:,:,:,:,:,1) = reshape(E,pts,pts,pts,pts,pts);
L2Score = (1-(L2/max(L2)));
quality(:,:,:,:,:,1) = reshape(L2Score,pts,pts,pts,pts,pts);

for i = 1:length(mag(:))
    if mag(i) < (104/320)
        magdummy(i) = 0;
    else
        magdummy(i) = ((1-104/320)^-1)*mag(i)-(((1-104/320)^-1)-1);
    end
end

quality(:,:,:,:,:,2) = reshape(magdummy,pts,pts,pts,pts,pts);
clear magdummy;

quality(:,:,:,:,:,3) = (0.5*tanh(10*(effRad./r1 -.96)*2*pi) + 0.5);
% quality(:,:,:,:,:,3) = 0.6107*(tanh(3*(effRad./r1 - .96)*2*pi) + 1);

n = size(mag,1);


q1 = sum(quality,6);
q1 = squeeze(q1);

q2 = q1(:);

v1 = zeros(length(mag)^5,5);
counter = 1;

var1 = [Lset; bEllSet; dlaserSet; eccSet; hset;]';

for i = 1:length(hset)
    for j = 1:length(eccSet)
        for k = 1:length(dlaserSet)
            for l = 1:length(bEllSet)
                for m = 1:length(Lset)
                    
                    v1(counter,1) = Lset(m);
                    v1(counter,2) = bEllSet(l);
                    v1(counter,3) = dlaserSet(j);
                    v1(counter,4) = eccSet(i);
                    v1(counter,5) = hset(k);
                    counter = counter +1;
                end
            end
        end
    end
end


stableq2 = zeros(size(q2));
y = 1; %how many pixels to step away
for z = 1:numSims
    
    [o1,o2,o3,o4,o5] = ind2sub([length(Lset),length(Lset),length(Lset),length(Lset),length(Lset)],z);
    
    Ci = [o1,o2,o3,o4,o5] - y;
    Cf = [o1,o2,o3,o4,o5] + y;
    
    if any(([o1,o2,o3,o4,o5] - y) < 1)
        
%         Ci = [o1,o2,o3,o4,o5];
        Ci([o1,o2,o3,o4,o5]-y < 1) = 1;
    end
        
    if any(([o1,o2,o3,o4,o5] + y) > length(Lset))
            
        Cf = [o1,o2,o3,o4,o5] + y;
        Cf(Cf>length(Lset)) = length(Lset);
    end
        
    
    
    set1 = zeros((1+2*y)^5,1); 
    counter = 1;
    for i = Ci(5):Cf(5)
        for j = Ci(4):Cf(4)
            for k = Ci(3):Cf(3)
                for l = Ci(2):Cf(2)
                    for m = Ci(1):Cf(1)

                        set1(counter) = q1(m,l,k,j,i);

                        counter = counter + 1;
                    end
                end
            end
        end
    end
    set1(set1 == 0) = [];
    stableq2(z) = min(set1);
    
end

figure(1)
binranges=[0:0.3:3];
a=histcounts(q1,[binranges Inf]);
b=histcounts(stableq2,[binranges Inf]);
bar(binranges,[a;b]')
xlabel('scores','FontSize',16,'FontWeight','bold')
xticks(0:0.3:3)
ylabel('counts','FontSize',16,'FontWeight','bold')
title('zoomed data set','FontSize',16)
legend('raw scores','stablized scores','Location','northeast','FontSize',12)

figure(2)
binranges=[2.82:0.005:2.89];
a=histcounts(q1,[binranges Inf]);
b=histcounts(stableq2,[binranges Inf]);
bar(binranges,[a;b]')
xlabel('scores','FontSize',16,'FontWeight','bold')
xticks(2.81:0.01:2.89)
ylabel('counts','FontSize',16,'FontWeight','bold')
title('zoomed data set','FontSize',16)
legend('raw scores','stablized scores','Location','northeast','FontSize',12)

% figure(2)
% binranges=[0:0.1:1];
% a=histcounts(quality(:,:,:,:,:,3),[binranges Inf]);
% b=histcounts(quality(:,:,:,:,:,2),[binranges Inf]);
% c=histcounts(quality(:,:,:,:,:,1),[binranges Inf]);
% bar(binranges,[a;b;c]')
% xlabel('scores','FontSize',16,'FontWeight','bold')
% xticks(binranges)
% ylabel('counts','FontSize',16,'FontWeight','bold')
% title('zoomed data set','FontSize',16)
% legend('effRad','mag','RMSE','Location','northwest','FontSize',12)