function [score, score_01, score_02, score_03] = Solid_Tank_Sim(x, gel_option)

%Variables
Lset = x(1);
hset = x(2);
dlaserSet = x(3);
bEllSet = x(4);
eccSet = x(5);


%SimSettings

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

N = 512;
r = 400000;
bath = 1.4886; %acrylic
span1 = 30*pi/180;
span2 = -30*pi/180;
scale = 1;
lensType = 2; %1 = circle, 2 = ellipse
radFaceSet = 50;
polAng = (90*pi/180);
%%% Attenuators
anglesSet = 0; %0:1:360 step size of rotation
theta = 0; %location of center point radially
radiusToImage = 0; %radius from center of bore to center of attenuator
r4 = 47.6250*scale; %radius of attenuator
diaset = 101.6;
d1set = 145+(diaset/2); % distance from front edge of acrylic block to center of bore
d2set = [275:25:400]-d1set; % distance from center of bore hole to detectors - 344mm = total length of acrylic block

%Simulation
a=1;
anglesSet = 0;

h = hset; %xaxis
k = 0;  %yaxis
L = Lset; %block length
radFace = radFaceSet(a);

bEll = bEllSet*scale;
dlaser = dlaserSet;
ecc = eccSet;
dia = diaset;
d1 = (L/2)+h;
d2 = L/2-h;


imangles=anglesSet;%find angle after each rotation
imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord
realh = 0; %this is so that the bore center appears at 0,0 graphically

[xInts,yInts,iRay,aRay] = RefracFunc_Merge(N,r,gel,span1,span2,scale,...
    dlaser,bath,dia,d1,d2,lensType,ecc,bEll,realh,k,radFace,polAng,...
    imCoordx,imCoordy,imangles,theta,r4, u);

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
        elseif xInts(i,j) == d2*scale
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
    if max(yInts(i,:))<((128+k)*scale) && min(yInts(i,:))>-((128+k)*scale) && nnz(xInts(i,:))>8 && max(xInts(i,:)) == d2*scale
        IntPts = [xInts(i,5), yInts(i,5); xInts(i,4), yInts(i,4)];
        shortRay(i,1) = pdist(IntPts,'euclidean');
    else
        shortRay(i,1) = 0;
    end
end
oopsies = (length(shortRay) - length(shortRay(shortRay>0)))/2;
[dumdum, effRadSelect] = min(shortRay(shortRay>0));
for i = effRadSelect + oopsies
    
    chordx = [xInts(i,5) xInts(i,4)];
    chordy = [yInts(i,5) yInts(i,4)];
    
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
for i = 1:320
    if proGeoPol(i) > max(proGeoPol)*.1
        magCount = magCount + 1;
        
    end
end

mag = (magCount/(r3*2))/(320/(r3*2));
%if mag < (104/320)
%    magdummy = 0;
%else
magdummy = ((1-104/320)^-1)*mag-(((1-104/320)^-1)-1);
%end

score_01 = -0.5*tanh(0.2*(beamUni - 5)*2*pi) + 0.5;
score_02 = magdummy;
score_03 = 0.5*tanh(10*(effRad/r1 -.96)*2*pi) + 0.5;
score = score_01 + score_02 + score_03;


end