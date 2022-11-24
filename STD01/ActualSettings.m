clear 
close all

N = 512;
r = 200;
ra1 = 400000;
pts = 1;

% gel = 1.4225;
gel = 1.34468 %Modus clearview(??) estimate (12% PG and 88% water) PG at 1.4384
% gel = 1.3319;
% gel = 1.4886;
bath = 1.4886; %acrylic
% bath = 1; %mismatch
% bath = 1.3319; %water
span1 = 30*pi/180;
span2 = -30*pi/180;


scale = 1;

% dlaserSet = 41; %actual
dlaserSet = 58; %fricke
% dlaserSet = 50 %flexy
% dlaserSet = 49; %w/ water
lensType = 2; %1 = circle, 2 = ellipse

% eccSet = 0; %actual
eccSet = 0; %fricke
% eccSet = 0; %flexy
% eccSet = 0.5;

% bEllSet = 61; %actual
bEllSet = 69; %fricke
% bEllSet = 66; %flexy
% bEllSet = 67;

radFaceSet = 0;
polAng = (90*pi/180);

%%% Attenuators
anglesSet = 0; %0:1:360 step size of rotation
theta = 0; %location of center point radially
radiusToImage = 0; %radius from center of bore to center of attenuator
r4 = 47.6250*scale; %radius of attenuator

%%% physical properties
diaset = 101.6; % 101.6mm - diameter of container
% hset = 4; %actual
hset = 13; %fricke
% hset = 6.5; %flexy
% hset = -33;

% Lset = 290; %actual
Lset = 204; %fricke
% Lset = 314; %flexy
% Lset = 448;