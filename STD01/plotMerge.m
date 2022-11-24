%Input the number of rays
clear
close all
ActualSettings
rIntensity = zeros(ra1,1);
a = 1;
a1 = 1;
a2 = 1;

% for z = 1:length(Lset)
for z = 1:1
    
    clf
    L = Lset(a);
    h = hset(a); %xaxis
    k = 0;  %yaxis
    radFace = radFaceSet(a)*scale;
    dlaser = dlaserSet(a);
    dia = diaset(a);
    d1 = (L/2)+h;
    d2 = L/2-h;
    ecc = eccSet(a);
    bEll = bEllSet(a);
    
    for boop = 1:length(anglesSet)
        imangles=anglesSet(boop);%find angle after each rotation
        imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
        imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord
        
        realh = 0; %this is so that the bore center appears at 0,0 graphically
        
        [xInts,yInts,iRay,aRay] = RefracFunc_Merge(N,r,gel,span1,span2,scale,...
        dlaser,bath,dia,d1,d2,lensType,ecc,bEll,realh,k,radFace,polAng,...
        imCoordx,imCoordy,imangles,theta,r4);
        hold on
        tic
        %Measured Values
        wallT = 0.125*25.4; %0.125in wall - container wall thickness
        gapT = (104 - dia)/2; %gap thickness 104 is bore hole diameter mm
        doff = 0;
        h = h*scale; k = k*scale;
        % dlaser = 5;
        
        %counting/conditional variables for loops
        rCount = 0;
        gCount = 0;
        bCount = 0;
        cCount = 0;
        gTrue = 0;  %start as false
        bTrue = 0;
        rTrue = 0;
        cTrue = 0;
        
        r1 = dia*scale/2 - wallT*scale; %gel radius
        r2 = dia*scale/2; %container radius
        r3 = dia*scale/2 + gapT*scale; %bore radius
        for i = 1:r
            for j = 1:8
                if xInts(i,j) ~= 0
                    x = [xInts(i,j) xInts(i,j + 1)];
                    y = [yInts(i,j) yInts(i,j + 1)];
                    if max(yInts(i,:))>(128*scale) || min(yInts(i,:))<-(128*scale)
                        l1 = line(x,y);
                        l1.Color = 'red';
                        rTrue = 1;
                        
                    elseif(nnz(xInts(i,:))<8 || max(xInts(i,:)) < d2*scale)
                        l2 = line(x,y);
                        if max(xInts(i,:)) < d2*scale
                            l2.Color = 'red';
                            rTrue = 1;
                        elseif nnz(xInts(i,:))<4
                            l2.Color = 'blue';
                            bTrue = 1;
                        else
                            l2.Color = 'cyan';
                            cTrue = 1;
                        end
                    else
                        l3 = line(x,y);
                        l3.Color = '[0,0.5,0]';
                        gTrue = 1;
                    end
                end
            end
            
            if rTrue == 0
                rIntensity(i) = iRay(i,find(iRay(i,:) > 0, 1));
            end
            
            if gTrue == 1
                gCount = gCount + 1;
                gTrue = 0;
            elseif bTrue == 1
                bCount = bCount + 1;
                bTrue = 0;
            elseif rTrue == 1
                rCount = rCount + 1;
                rTrue = 0;
            elseif cTrue == 1
                cCount = cCount + 1;
                cTrue = 0;
            end
        end
        
        theta = 0:0.001:2*pi;        
    
        %calc chord for effective rad
        %first green ray
        for i = 1:r
            if max(yInts(i,:))<(128*scale) && min(yInts(i,:))>-(128*scale) && nnz(xInts(i,:))>8 && max(xInts(i,:)) == d2*scale
                chordx = [xInts(i,5) xInts(i,4)];
                chordy = [yInts(i,5) yInts(i,4)];
                break
            end
        end
        
        %draw over first green ray
        line(chordx,chordy,'LineWidth',2,'Color','m')
        %math for finding crossing point between center and line segment
        slope = (diff(chordx)./diff(chordy))^-1;
        yIntercept = chordy(1)-slope*chordx(1); %yint of chord
        yIntercept2 = k-(-slope^-1)*realh; %yint of line from center
        xIntersect = (yIntercept2-yIntercept)/(slope - (-slope^-1));
        yIntersect = slope*(xIntersect)+yIntercept;
        %plot intersect point
        plot(xIntersect,yIntersect,'black*','markersize',8)
        %plot line from center to intersect point
        line([realh,xIntersect],[k,yIntersect],'LineWidth',2)
        %calc effRad
        effRad = norm([xIntersect-realh,yIntersect-k])
        
        % line([-d1*scale d2*scale],[0 0]);
        %     axis([-1000 500 -500 500])
        % axis equal
        title ('Ray Plot - acrylic - acrylic (no dye)');
        %     lgd = legend([l1,l3,l2],'Undetected Rays','Attenuated Ray','Direct-path Rays');
        lgd.Location = 'northwest';
        axis square
        

    end
end

x1 = -d1*scale;
if lensType == 1 %circle
    hFace = (2*x1 + ((2*x1)^2 - 4*(x1^2 - radFace^2))^0.5)/2;
    kFace = 0;
    %circle
    plot(hFace+radFace*sin(theta(3142:6284)),kFace+radFace*cos(theta(3142:6284)),'LineWidth' ,1, 'Color', 'blue');
    
elseif lensType == 2 %ellipse
    hFace = x1+(bEll^2/(1 + ecc^2))^0.5; %b term from ellipse equation using a and ecc
    kFace = 0;
    %ellipse
    plot(hFace+(bEll^2/(1 + ecc^2))^0.5*cos(theta(1571:4713)),kFace+bEll*sin(theta(1571:4713)),'LineWidth' ,1, 'Color', 'cyan');
    
end

plot(realh+r1*sin(theta),k+r1*cos(theta), 'LineWidth',1, 'Color', 'blue');
plot(realh+r2*sin(theta),k+r2*cos(theta), 'LineWidth',1, 'Color', 'm');
plot(realh+r3*sin(theta),k+r3*cos(theta), 'LineWidth',1, 'Color', 'black');
line([-d1*scale -d1*scale],[150*scale -150*scale]);
line([d2*scale d2*scale],[128*scale -128*scale],'LineWidth',4,'Color','[0.85 0.325 0.098]');
line([d2*scale d2*scale],[150*scale -150*scale]);
line([min(min(xInts))*scale min(min(xInts))*scale],[10*scale -10*scale]);
line([min(min(xInts))*scale (min(min(xInts))-10)*scale],[10*scale 10*scale]);
line([min(min(xInts))*scale (min(min(xInts))-10)*scale],[-10*scale -10*scale]);
xlim([-300 300])
ylim([-300 300])
hold off

for z = 1:1

    L = Lset(a);
    h = hset(a); %xaxis
    k = 0;  %yaxis
    radFace = radFaceSet(a)*scale;
    dlaser = dlaserSet(a);
    dia = diaset(a);
    d1 = (L/2)+h;
    d2 = L/2-h;
    ecc = eccSet(a);
    bEll = bEllSet(a);
    
    for boop = 1:length(anglesSet)
        imangles=anglesSet(boop);%find angle after each rotation
        imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
        imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord
        
        realh = 0; %this is so that the bore center appears at 0,0 graphically
        
        [xInts,yInts,iRay,aRay] = RefracFunc_Merge(N,ra1,gel,span1,span2,scale,...
        dlaser,bath,dia,d1,d2,lensType,ecc,bEll,realh,k,radFace,polAng,...
        imCoordx,imCoordy,imangles,theta,r4);
        hold on
        tic
    end
    
        %finding first none zero intesity values for each row the reachs detector
    %ie the value at detector
    I = zeros(ra1,1);
    Y = zeros(ra1,1);
    for i = 1:ra1
        if nnz(xInts(i,:)) < 8 %filter out tangent rays for profile
        else
            for j = 1:8
                if (xInts(i,j) == d2*scale)
                    I(i,1)=iRay(i,j+1);
                    Y(i,1)=yInts(i,j);
                end
            end
        end
    end

    
end
    
proGeo = zeros(1,320);
proGeoPol = zeros(1,320);
for i = 1:ra1
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


figure(2)
% plot(proGeo)
hold on
plot(proGeoPol,'m','LineWidth', 3)
title ({'400,000 ray SimProfile';'Local Maxima Beam Uniformity'},'FontSize',14,'FontWeight','bold');
axis square
yline(mean(proGeoPol),'--b','LineWidth', 3)
ylim([0 5000])
% ylim([0 2000])
% ylim([0 800])
xlim([0 320])
xlabel('Detector #','FontSize',14,'FontWeight','bold')
ylabel('Intensity','FontSize',14,'FontWeight','bold')
legend('simulated profile','mean','Location','north','FontSize',12)

