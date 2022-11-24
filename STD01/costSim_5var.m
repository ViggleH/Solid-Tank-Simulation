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

effRad = zeros(1,length(bEllSet));
a=1;
counter = 0;

anglesSet = 0;

for v1 = 1:length(Lset)
    for v2 = 1:length(bEllSet)
        for v3 = 1:length(hset)
            for v4 = 1:length(dlaserSet)
                for v5 = 1:length(eccSet)
                    
                    tS =tic;
                    
                    h = hset(v3); %xaxis
                    k = 0;  %yaxis
                    L = Lset(v1); %block length
                    radFace = radFaceSet(a);
                    
                    bEll = bEllSet(v2)*scale;
                    dlaser = dlaserSet(v4);
                    ecc = eccSet(v5);
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
                    tanRays = zeros(r,1);
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
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %plotting dynamic range profiles
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %             figure(2)
                    %             %     plot(proGeo)
                    %             title('Dynamic Range (Profile)')
                    %             hold on
                    %             plot(proGeoPol)
                    %             legend(sprintfc('circle Radius %d mm',bEllSet),'Location','westoutside')
                    noZpro = nonzeros(proGeoPol');
                    beamUni(v1,v2,v3,v4,v5) = max(proGeoPol)/prctile(noZpro,10);
                    
                    
                    %         figure(3)
                    %         title('Dynamic Range (max/min)')
                    %         hold on
                    %         bar(bEllSet(z),max(proGeoPol)/min(proGeoPol))
                    %
                    %             figure(4)
                    %             title('Profile Ratio (Polarity/Geometry)')
                    %             plot(proGeoPol./proGeo)
                    %             hold on
                    
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
                    effRad(v1,v2,v3,v4,v5) = norm([xIntersect-realh,yIntersect-k]);
                    
                    %magnification, how many detectors are viewing the r3 radius
                    %ideal mag is 104mm viewed by 320 detectors = 3.08x
                    magCount = 0;
                    for i = 1:320
                        if proGeoPol(i) > max(proGeoPol)*.1
                            magCount = magCount + 1;
                            
                        end
                    end
                    
                    mag(v1,v2,v3,v4,v5) = (magCount/(r3*2))/(320/(r3*2));
                    
                    counter = counter +1;
                    sprintf('Done with set %d of %d after %3.2f seconds',counter, length(Lset)*length(bEllSet)*length(hset)*length(dlaserSet)*length(eccSet), toc(tS))
                end
            end
        end
    end
end
toc(tStart)