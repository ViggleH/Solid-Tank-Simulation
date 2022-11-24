%##########################################################################
%                       OPTIMAL SCANNER GEOMETRY
%##########################################################################
%All geometry values are in mm and not scaled

function [geo] = ScannerGeometry_Optimal(r,numDet,laserType,N)


geo.N = N;

geo.r = r;
geo.numDet = numDet;

geo.laserType = laserType;

geo.arrayW = 35;

geo.dlaser = 41; %optimal
% geo.dlaser = 41; %from systemAligner
% geo.dlaser = 50;
% geo.doff = 0; %optimal
% geo.doff = 1.75; %from systemAligner
% geo.doff = 1.3; %2021 alignment
% geo.doff = -1.15; %spatial res alignment
% geo.doff = 0.95; %spatial res
geo.doff = 0;
% geo.doff = -0.35; %needle alignment
% geo.doff = 1;
% geo.doff = 0.35
% geo.hh = 0;

% geo.h = 4; %optimal
geo.h = 40;
geo.k = 0; %optimal
% geo.k = -7.0;

% geo.bath = 1.4886;
geo.bath = 1.3319;
geo.gel = 1.4225; %optimal
% geo.gel = 1.3319; %water
% geo.gel = 1.00;
geo.air = 1.00029; % refractive index of air
geo.acrylic = 1.4886; % refractive index of acrylic (block and container)

geo.span1 = 30*pi/180; %optimal
% geo.span1 = 23*pi/180; %from systemAligner
geo.span2 = -30*pi/180; %optimal
% geo.span2 = -26*pi/180; %from systemAligner

geo.dia = 101.6; % diameter of container
geo.wallT = 0.125*25.4; % %0.125in - container wall thickness
geo.gapT = (104-geo.dia)/2; % gap thickness i.e. bath thickness

geo.r1 = (geo.dia/2) - geo.wallT; %gel radius
geo.r2 = (geo.dia/2); %container radius
geo.r3 = (geo.dia/2) + geo.gapT; %bore radius

geo.lensType = 'ellipse'; % Block face lens type (circle or ellipse)
geo.radFace = 61; % Block face radius (if circle)
geo.bEll = 61; % Block face semi-major axis length (if ellipse)
geo.ecc = 0; % Block face eccentricity (if ellipse)
geo.L = 290; % Total length of block
% geo.L = 250;

geo.bEll2 = 65; %semi-major axis rear wall
geo.ecc2 = 1; % Ecc of the rear wall

geo.polAng = (90*pi/180);

geo.d1 = (geo.L/2) + geo.h; % distance from front edge of acrylic block to center of bore
geo.d2 = (geo.L/2) - geo.h; % distance from center of bore hole to detectors, L = total length of acrylic block

geo.det = geo.d2+200; % detectors position

geo.detBayHeight = 320*0.8;

%% Andy Edit - h and k correspond to the bore center position and influence x0 and y0 because of a 
% desire to keep bore center at 0,0
% switch(geo.laserType)
%   case('optimalCrossover')
%     geo.x0 = -geo.d1 - geo.dlaser + geo.h; % laser source position on x-axis 
%     geo.y0 = geo.doff + geo.k; % laser source position on y-axis
%   case('fan')
%     geo.x0 = -geo.d1 - geo.dlaser + geo.h; % laser source position on x-axis 
%     geo.y0 = geo.doff + geo.k; % laser source position on y-axis
%   case('parallel')
%     geo.x0 = ones(1,geo.r) * (-geo.d1 - geo.dlaser + geo.h);
%     geo.y0 = linspace((geo.doff + geo.k) - (geo.dia/2) , (geo.doff + geo.k) + (geo.dia/2) ,r);
%   case('crossover')
%     geo.x0 = ones(1,geo.r) * (-geo.d1 - geo.dlaser + geo.h);
%     geo.y0 = linspace((geo.doff + geo.k) - (geo.dia/2) , (geo.doff + geo.k) + (geo.dia/2) ,geo.r);
%     
%     geo.bottomDetYcross = 1;    
% end

%Andy Edit
switch(geo.laserType)
  case('optimalCrossover')
    geo.x0 = -geo.d1 - geo.dlaser; % laser source position on x-axis 
    geo.y0 = geo.doff; % laser source position on y-axis
  case('fan')
    geo.x0 = -geo.d1 - geo.dlaser; % laser source position on x-axis 
    geo.y0 = geo.doff; % laser source position on y-axis
  case('parallel')
    geo.x0 = ones(1,geo.r) * (-geo.d1 - geo.dlaser);
    geo.y0 = linspace((geo.doff) - (geo.dia/2) , (geo.doff) + (geo.dia/2) ,r);
  case('crossover')
    geo.x0 = ones(1,geo.r) * (-geo.d1 - geo.dlaser);
    geo.y0 = linspace((geo.doff) - (geo.dia/2) , (geo.doff) + (geo.dia/2) ,geo.r);
    
    geo.bottomDetYcross = 1;    
end

geo.wall = -geo.d1; % position of front edge of acrylic block

assignin('base', 'r1', geo.r1)
assignin('base', 'r2', geo.r2)
assignin('base', 'r3', geo.r3)
 
return

end



