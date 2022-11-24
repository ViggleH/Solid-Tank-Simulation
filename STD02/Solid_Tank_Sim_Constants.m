%SimSettings
N = 512;
r = 400000;
% r = 80;

arrayW = 50;

%gel = 1.4225; %FlexyDos3D
% gel = 1.34468; %Modus clearview(??) estimate (12% PG and 88% water) PG at 1.4384

% bath = 1.4886; %acrylic
% bath = 1.3319; %water
bath = 1.467; %parafin oil
air = 1.00029; % refractive index of air
acrylic = 1.4886; % refractive index of acrylic (block and container)

span1 = 30*pi/180;
span2 = -30*pi/180;
doff = 0;
k = 0;

scale = 1;

lensType = 2; %1 = circle, 2 = ellipse

radFaceSet = 50;
polAng = (90*pi/180);
%%% Attenuators

anglesSet = 0; %0:1:360 step size of rotation
theta = 0; %location of center point radially
radiusToImage = 0; %radius from center of bore to center of attenuator
imangles=anglesSet;%find angle after each rotation

r4 = 47.6250*scale; %radius of attenuator

dia = 101.6; % diameter of container
wallT = 0.125*25.4; % %0.125in - container wall thickness
gapT = (104-dia)/2; % gap thickness i.e. bath thickness

r1 = (dia/2) - wallT; %gel radius
r2 = (dia/2); %container radius
r3 = (dia/2) + gapT; %bore radius
lensType = 'ellipse';
laserType = 'line';
