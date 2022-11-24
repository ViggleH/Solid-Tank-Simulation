%##########################################################################
%                     Intesity plot
%##########################################################################
% this program plots the intesity of the light rays as they pass through
% the block and stuff, with a histogram of the
% intesity being the number of counts with the Yints being the x-axis
clear
close all
SimSettings
tStart = tic;

effRad = zeros(1,length(radFaceSet));
a=1;
counter = 0;

anglesSet = 0;

var1 = Lset;
var2 = bEllSet;
var3 = dlaserSet;
var4 = eccSet;
var5 = hset;

% simSpace is the size of the space you're exploring
simSpace = [length(var1), length(var2), length(var3), length(var4), length(var5)];

parfor_progress(prod(simSpace));

% calculate the total number of simulations
numSims  = prod(simSpace);
% pre-allocate data
noZbeamUni10 = zeros(numSims,1);
noZbeamUni5 = zeros(numSims,1);
effRad = zeros(numSims,1);
mag = zeros(numSims,1);
data = zeros(numSims,1);
profile = cell(numSims,1);


parfor idx = 1:numSims
    % convert scalar index into subscripts
    [v1, v2, v3, v4, v5] = ind2sub(simSpace, idx);
    % Index input vector (you should be able to use a single vector
    % subscript here, as shown)
    
    
    
    h = hset(v5); %xaxis
    k = 0;  %yaxis
    L = Lset(v1); %block length
    radFace = radFaceSet(a);
    
    bEll = bEllSet(v2)*scale;
    dlaser = dlaserSet(v3);
    ecc = eccSet(v4);
    dia = diaset;
    d1 = (L/2)+h;
    d2 = L/2-h;
    
    imangles=anglesSet;%find angle after each rotation
    imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
    imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord
    realh = 0; %this is so that the bore center appears at 0,0 graphically
    
    [xInts,yInts,iRay,aRay] = RefracFunc_Merge(N,r,gel,span1,span2,scale,...
        dlaser,bath,dia,d1,d2,lensType,ecc,bEll,realh,k,radFace,polAng,...
        imCoordx,imCoordy,imangles,theta,r4);
    
    h = h*scale; k = k*scale;
    
    
        det = d2*scale;
    wallT = 0.125*25.4; %0.125in wall - container wall thickness
    gapT = (104 - dia)/2; %gap thickness 104 is bore hole diameter mm
    doff = 0; % vertical offset of laser source
    
    r1 = diaset*scale/2 - wallT*scale; %gel radius
    r2 = diaset*scale/2; %container radius
    r3 = diaset*scale/2 + gapT*scale; %bore radius
    
    
    %finding first none zero intesity values for each row the reachs detector
    %ie the value at detector
    I = zeros(r,1);
    Y = zeros(r,1);
    for i = 1:r
        if nnz(xInts(i,:)) < 8 %filter out tangent rays for profile
        else
            for j = 1:8
                if (xInts(i,j)==det )
                    I(i,1)=iRay(i,j+1);
                    Y(i,1)=yInts(i,j);
                end
            end
        end
    end
    
    
    proGeo = zeros(1,320);
    proGeoPol = zeros(1,320);
    for i = 1:r
        for j = 1:8
            if nnz(xInts(i,:)) < 8 %filter out tangent rays for profile
            elseif xInts(i,j) == d2*scale;
                for f = 1:320
                    bin1 = (128 - (f-1)*0.8)*scale;
                    bin2 = (128 - f*0.8)*scale;
                    if  (bin1 >= yInts(i,j)) && (yInts(i,j) >= bin2)
                        proGeo(1,f) = proGeo(1,f) + 1;
                        proGeoPol(1,f) = proGeoPol(1,f) + I(i);
                    end
                end
            end
        end
    end
    
    
    
    %     dumVar = sort(proGeoPol(1:length(proGeoPol)/2));
    %     beamUni(idx) = dumVar(length(proGeoPol)/2*.9)/dumVar(length(proGeoPol)/2*.1);
    %     clear dumVar
    
    profile{idx} = proGeoPol;
    
    noZpro = nonzeros(proGeoPol');
    %     noZpro = reshape(goop,pts,pts,pts,pts)'
    
    %     beamUni10(idx) = max(proGeoPol)/prctile(proGeoPol,10);
    %     beanUni5(idx) = max(proGeoPol)/prctile(proGeoPol,5);
    
    noZbeamUni10(idx) = max(proGeoPol)/prctile(noZpro,10);
    noZbeamUni5(idx) = max(proGeoPol)/prctile(noZpro,5);
    
    counts = histcounts(proGeoPol,[1:1000]*100);
    p = counts/sum(counts);
    U(idx) = sum(p.^2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calc cord for effective rad
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %first green ray
    for i = 1:r
        if max(yInts(i,:))<(128*scale) && min(yInts(i,:))>-(128*scale) && nnz(xInts(i,:))>8 && max(xInts(i,:)) == d2*scale
            chordx = [xInts(i,5) xInts(i,4)];
            chordy = [yInts(i,5) yInts(i,4)];
            break
        end
    end
    %math for finding crossing point between center and line segment
    slope = (diff(chordx)./diff(chordy))^-1;
    yIntercept = chordy(1)-slope*chordx(1); %yint of chord
    yIntercept2 = k-(-slope^-1)*realh; %yint of line from center
    xIntersect = (yIntercept2-yIntercept)/(slope - (-slope^-1));
    yIntersect = slope*(xIntersect)+yIntercept;
    
    %calc effRad
    effRad(idx) = norm([xIntersect-realh,yIntersect-k]);
    
    %magnification, how many detectors are viewing the r3 radius
    %ideal mag is 104mm viewed by 320 detectors = 3.08x
    magCount = 0;
    for i = 1:320
        if proGeoPol(i) > mean(proGeoPol)*.1
            magCount = magCount + 1;
            
        end
    end
    
    mag(idx) = (magCount)/(320);
    
    data(idx) = v1;


    parfor_progress;
end

noZbeamUni10 = reshape(noZbeamUni10,pts,pts,pts,pts,pts);
noZbeamUni5 = reshape(noZbeamUni5,pts,pts,pts,pts,pts);
effRad = reshape(effRad,pts,pts,pts,pts,pts);
mag = reshape(mag,pts,pts,pts,pts,pts);
data = reshape(data,pts,pts,pts,pts,pts);
U = reshape(U,pts,pts,pts,pts,pts);

parfor_progress(0);

% save test.mat
save june21_fricke8.mat

toc(tStart)
    
    
    
    
    
    
    
    
    
    
    
    