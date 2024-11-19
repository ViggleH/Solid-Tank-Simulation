%##########################################################################
%                       OPTIMAL SCANNER GEOMETRY
%##########################################################################
%All geometry values are in mm and not scaled

function [geo] = setup_geometry(r,numDet,laserType,N, x)


geo.N = N;

geo.r = r;
geo.numDet = numDet;

geo.laserType = laserType;

% geo.arrayW = 30;
geo.arrayW = 15;

% geo.dlaser = 41; %optimal
% geo.dlaser = 41; %from systemAligner
geo.dlaser = 100;
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
geo.h = 0;
geo.k = 0; %optimal
% geo.k = -7.0;

% geo.bath = 1.4886; %acrylic matching fluid
% geo.bath = 1.3319; %water
geo.bath = 1.467; %parafin oil
% geo.bath = 1.00029; %air
% geo.gel = 1.4225; %optimal
geo.gel = 1.4225;
% geo.gel = 1.34468; %clearview
% geo.gel = 1.3319; %water
% geo.gel = 1.00;
geo.air = 1.00029; % refractive index of air
geo.acrylic = 1.4886; % refractive index of acrylic (block and container)
% geo.block = geo.air;
geo.block = geo.acrylic;

geo.span1 = 16.35*pi/180; %optimal
% geo.span1 = 12*pi/180; %from systemAligner
geo.span2 = -16.35*pi/180; %optimal
% geo.span2 = -12*pi/180; %from systemAligner

% geo.dia = 152.4; % 6"
% geo.dia = 127; % 5"
geo.dia = 101.6; % 4" diameter of container
% geo.dia = 3.5*25.4; % 3.5"
% geo.dia = 76.2; % 3"
% geo.dia = 50.8; % 2"

geo.wallT = 0.125*25.4; % 0.125in - container wall thickness
% geo.wallT = 0.25*25.4; % 0.25in
% geo.wallT = 0.5*25.4; % 0.5in


% geo.gapT = (104-geo.dia)/2; % gap thickness i.e. bath thickness
geo.gapT = 1.2;

geo.r1 = (geo.dia/2) - geo.wallT; %gel radius
geo.r2 = (geo.dia/2); %container radius
geo.r3 = (geo.dia/2) + geo.gapT; %bore radius

geo.lensType = 'ellipse'; % Block face lens type (circle or ellipse)
geo.radFace = 61; % Block face radius (if circle)
% geo.bEll = 61; % Block face semi-major axis length (if ellipse)
geo.bEll = 160;
% geo.ecc = 0; % Block face eccentricity (if ellipse)
geo.ecc = 1; 
% geo.L = 290; % Total length of block
geo.L = 114;

% geo.bEll2 = 65; %semi-major axis rear wall
geo.bEll2 = 50; 
% geo.ecc2 = 1; % Ecc of the rear wall
geo.ecc2 = 0; % Ecc of the rear wall

geo.polAng = (90*pi/180);
% d3 = 18.6; %distance from back ellipse to detectors


%dom check
% x = [200.5 -14.8 70.0 91.3 1.0 70 0 130.0] %andy guess
% x = [114 0 100 67.8 0 82.1 0 10.3] %dom value
%x = [114 0 100 67.8 0 80 0 200]; %andy improve
% x =   [114 0 40 148.62 0.6754 93.4234 2.2803 286]
% x = [110 0 40 100 0 100 0 60] %andy air scan tests

%x = [128.7352 3.7393 150 84.975 0 73.3395 0 88];
% x = [140.7352 -2.2607 150 84.975 0 73.3395 0 88]; %thicker
% x = [140.7352 -2.2607 150 84.975 0 73.3395 0 188];
% x = [140.7352 -2.2607 150 84.975 0 73.3395 0 118]; %thicker
% x =  [194 -40 100 77.0708 1 91.8702 0.4343 89.9055] %arrayW = 14
% x =  [124 0 100 77.0708 1 63.8702 0.4343 89.9055] %arrayW = 14

geo.L = x(1);
geo.h = x(2);
geo.dlaser = x(3);
geo.bEll = x(4);
geo.ecc = x(5);
geo.bEll2 = x(6);
geo.ecc2 = x(7);
geo.d3 = x(8);

geo.d1 = (geo.L/2) + geo.h; % distance from front edge of acrylic block to center of bore
geo.d2 = (geo.L/2) - geo.h; % distance from center of bore hole to detectors, L = total length of acrylic block


geo.det = geo.d2+geo.d3; % detectors position

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
    case('line')
        geo.x0 = -geo.d1 - geo.dlaser; % laser source position on x-axis
        geo.y0 = geo.doff; % laser source position on y-axis
end

geo.wall = -geo.d1; % position of front edge of acrylic block

assignin('base', 'r1', geo.r1)
assignin('base', 'r2', geo.r2)
assignin('base', 'r3', geo.r3)

 
return

end
