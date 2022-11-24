function [scoreSum, score1, score2, score3, score4] = Solid_Tank_Sim_Andy2022(x, gel_option)
%Constants
Solid_Tank_Sim_Constants

%Variables
L = x(1);
h = x(2);
dlaser = x(3);
bEll = x(4);
ecc = x(5);
bEll2 = x(6);
ecc2 = x(7);
d3 = x(8);

%Simulation
anglesSet = 0;

switch gel_option
    case 1 %water
        u = 0;
        gel = 1.3316;
    case 2 %Flexidos3d
        u = 0.0166;
        gel = 1.4225;
    case 3 %ClearView
        u = 0.0110;
        gel = 1.34468;
end


d1 = (L/2)+h;
d2 = L/2-h;

imangles=anglesSet;%find angle after each rotation
imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord
realh = 0; %this is so that the bore center appears at 0,0 graphically

[xInts,yInts,iRay,aRay] = RefracFunc_CCD_Sim(x, gel, u);

h = h*scale; k = k*scale;

det = d2+d3;
wallT = 0.125*25.4; %0.125in wall - container wall thickness
gapT = (104 - dia)/2; %gap thickness 104 is bore hole diameter mm
doff = 0; % vertical offset of laser source

r1 = dia/2 - wallT; %gel radius
r2 = dia/2; %container radius
r3 = dia/2 + gapT; %bore radius


%finding first none zero intesity values for each row the reachs detector
%ie the value at detector
I = zeros(r,1);
Y = zeros(r,1);
for i = 1:r
    if nnz(xInts(i,:)) < 9 %filter out tangent rays for profile
    else
        for j = 1:9
            if (xInts(i,j)==det )
                I(i,1)=iRay(i,j+1);
                Y(i,1)=yInts(i,j);
            end
        end
    end
end

detN = 200;
proGeo = zeros(1,detN);
proGeoPol = zeros(1,detN);
proNoCross = zeros(1,detN);
proCross = zeros(1,detN);
crossCount = 0;
for i = 1:r
    for j = 1:9
        if nnz(xInts(i,:)) < 9 %filter out tangent rays for profile
        elseif xInts(i,j) == det
            for f = 1:detN %overall profile
                bin1 = (arrayW - (f-1)*(2*arrayW/detN));
                bin2 = (arrayW - f*(2*arrayW/detN));
                if  (bin1 >= yInts(i,j)) && (yInts(i,j) >= bin2)
                    proGeo(1,f) = proGeo(1,f) + 1;
                    proGeoPol(1,f) = proGeoPol(1,f) + I(i);
                    if i <= r/2
                        if nnz(yInts(i,:))>8 && yInts(i,1)<yInts(i+1,1)
                            proCross(1,f) = proCross(1,f) + I(i);
                            crossCount = crossCount + 1;
                        else
                            proNoCross(1,f) = proNoCross(1,f) + I(i);
                        end
                    elseif i > r/2
                        if nnz(yInts(i,:))>8 && yInts(i,1)>yInts(i-1,1)
                            proCross(1,f) = proCross(1,f) + I(i);
                            crossCount = crossCount + 1;
                        else
                            proNoCross(1,f) = proNoCross(1,f) + I(i);
                        end
                    end
                end
            end
        end
    end
end

%figure(1)
%plot(proGeoPol)
%hold on
%plot(proNoCross)
%plot(proCross)
%legend('profile','primary','crossover')
%hold off

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
beamUni = max(proGeoPol)/prctile(noZpro,10);


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
%calc chord for effective rad
%first green ray
shortRay = zeros(r,1);
for i = 1:r
    if max(yInts(i,:))<(arrayW+k) && min(yInts(i,:))>-(arrayW+k) && nnz(xInts(i,:))>9 && max(xInts(i,:)) == det
        IntPts = [xInts(i,6), yInts(i,6); xInts(i,5), yInts(i,5)];
        shortRay(i,1) = pdist(IntPts,'euclidean');
    else
        shortRay(i,1) = 0;
    end
end
oopsies = (length(shortRay) - length(shortRay(shortRay>0)))/2;
[dumdum, effRadSelect] = min(shortRay(shortRay>0));
for i = effRadSelect + oopsies
    
    chordx = [xInts(i,6) xInts(i,5)];
    chordy = [yInts(i,6) yInts(i,5)];
    
end

%math for finding crossing point between center and line segment
if size(chordx, 1) > 0 && size(chordy, 1)
slope = (diff(chordx)./diff(chordy))^-1;
yIntercept = chordy(1)-slope*chordx(1); %yint of chord
yIntercept2 = k-(-slope^-1)*realh; %yint of line from center
xIntersect = (yIntercept2-yIntercept)/(slope - (-slope^-1));
yIntersect = slope*(xIntersect)+yIntercept;

%calc effRad
effRad = norm([xIntersect-realh,yIntersect-k]);
else
    effRad = 0;
end

%magnification, how many detectors are viewing the r3 radius
%ideal mag is 104mm viewed by 320 detectors = 3.08x
magCount = 0;
for i = 1:detN
    if proGeoPol(i) > max(proGeoPol)*.1
        magCount = magCount + 1;
        
    end
end

mag = magCount/detN;
%if mag < (104/320)
%    magdummy = 0;
%else
magdummy = ((1-104/320)^-1)*mag-(((1-104/320)^-1)-1);
%end

score1 = (-0.5*tanh(0.2*(beamUni - 5)*2*pi) + 0.5) * 0.5;
score2 = magdummy * 0.25;
score3 = 0.5*tanh((10/3)*(effRad/r1 -.88)*2*pi) + 0.5;
score4 = exp(5 * sum(proNoCross)/sum(proGeoPol) - 5) * 0.25;
scoreSum = score1+score2+score3+score4;


end