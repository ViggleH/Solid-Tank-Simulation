clear 
close all

N = 512;
r = 40000;
pts = 9;

gel = 1.4225; %Flexidos3d
%gel = 1.34468 %Modus clearview(??) estimate (12% PG and 88% water) PG at 1.4384
%gel = 1.3316; %Refraction index for water 635nm.


bath = 1.4886; %acrylic
% bath = 1.5; %mismatch
% bath = 1.3319; %water
span1 = 30*pi/180;
span2 = -30*pi/180;
% scale = 200/104*2;
scale = 1;
dlaserSet = linspace(52,56,pts);
% dlaserSet = linspace(35,74,pts);
% dlaserSet = [23 32 48 57];
% dlaserSet = 41;
lensType = 2; %1 = circle, 2 = ellipse
% eccSet = 0.5667; %eccentricity
% eccSet = linspace(0,0.9,pts);
eccSet = linspace(0.0,0.4,pts);
% bEllSet = 70.6667; %ellipse major axis (yaxis length)
% bEllSet = linspace(40,100,pts);
bEllSet = linspace(64,68,pts);
% radFaceSet = 40:7.5:100;
% radFaceSet = linspace(70,100,7);
% radFaceSet = [50 60 80 90];
radFaceSet = 50;
polAng = (90*pi/180);

%%% Attenuators
anglesSet = 0; %0:1:360 step size of rotation
theta = 0; %location of center point radially
radiusToImage = 0; %radius from center of bore to center of attenuator
r4 = 47.6250*scale; %radius of attenuator

%%% physical properties
diaset = 101.6; % 106mm - diameter of container
% hset = 8.3333; %bore translation mm
hset = linspace(7,11,pts);
% hset = [-27 -13 13 27];
% Lset = 320; %length of block mm
% Lset = linspace(300,380,pts);
Lset = linspace(201,209,pts);
% Lset = 325;


 d1set = 145+(diaset/2); % distance from front edge of acrylic block to center of bore
% % d1set = 75+(diaset/2);
 d2set = [275:25:400]-d1set; % distance from center of bore hole to detectors - 344mm = total length of acrylic block
% % d2set = 400-d1set;