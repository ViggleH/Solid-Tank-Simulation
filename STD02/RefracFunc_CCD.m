%##########################################################################
%          COMPLETE REFRACTION FUNCTION IS DETECTED ZERO DUPLICATES
%##########################################################################

function [XintersectionMatrix,YintersectionMatrix,rayNumber] = CompleteRefractionFunctionIsDetected_ZeroDuplicates(geo)
    
tStart = tic;
%################## Physical parameters of system #########################
%Can print system paramters for looping

% fprintf('r = %f\n' , r); % number of rays
% fprintf('gel = %f\n' , gel); % refractive index of the gel
% fprintf('span1 = %f\n' , span1); % the span of the fan beam above the horizontal(in radian) (Positive Entry)
% fprintf('span2 = %f\n' , span2); % the span of the fan beam below the horazontal(in radian) (Negative Entry)
% fprintf('edgeDensity = %f\n' , edgeDensity); % Fraction of total ray density incident on the edges of fan beam(.2 radian). 
%                                              % Uniform is edgeDesity = (.02+.02)/(span1+(-span2))
% fprintf('dlaser = %f\n' , dlaser); % distance from laser source to acrylic block
% fprintf('doff = %f\n' , doff); % vertical offset of laser source 
% fprintf('bath = %f\n' , bath); % refractive index of the matching bath
% fprintf('edgeSpan = %f\n' , edgeSpan); %angle(rad) of increased density on both ends of fan beam

r = geo.r; % number of rays
numDet = geo.numDet; % number of detectors

dlaser = geo.dlaser; % distance from laser source to acrylic block
doff = geo.doff; % vertical offset of laser source 

bath = geo.bath; % refractive index of the matching bath
gel = geo.gel; % refractive index of the gel
air = geo.air;% refractive index of air
acrylic = geo.acrylic; % refractive index of acrylic (block and container)

span1 = geo.span1; % the span of the fan beam above the horizontal(in radian) (Positive Entry)
span2 = geo.span2; % the span of the fan beam below the horazontal(in radian) (Negative Entry)

dia = geo.dia; % 104mm - diameter of bore hole
wallT = geo.wallT; % 3.5mm - container wall thickness
% hh = geo.hh;

h = geo.h; k = geo.k; % center of bore at (h,k)
h = 0;

gapT = geo.gapT; % gap thickness i.e. bath thickness
d1 = geo.d1; % distance from front edge of acrylic block to center of bore
d2 = geo.d2; % distance from center of bore hole to detectors, 344mm = total length of acrylic block
x0 = geo.x0; % laser source position on x-axis 
det = geo.det; % detectors position
y0 = geo.y0; % laser source position on y-axis
wall = geo.wall; % position of front edge of acrylic block

bEll2 = geo.bEll2;
ecc2 = geo.ecc2;

switch(geo.laserType)
  case 'fan'
    angles = linspace(span1,span2,r); % angle of inclination of each ray from laser source
  case 'optimalCrossover'
     angles = [linspace(-40*pi/180,span2,r/2) linspace(span1,40*pi/180,r/2)]; % angle of inclination of each ray from laser source
end

r1 = geo.r1; %gel radius
r2 = geo.r2; %container radius
r3 = geo.r3; %bore radius


% critical angles assuming bath < acrylic
crit1 = asin(bath/acrylic); % critical angle where total internal reflection occurs at Acrylic->Bath interface
crit2 = asin(gel/acrylic); % critical angle where total internal reflection occurs at Gel->Acrylic interface
crit3 = asin(air/acrylic); % critical angle where total internal reflection occurs at Acrylic->air interface

%Calculate position of middle of each detector
count = 1;
detSpace = geo.detBayHeight/geo.numDet; %Detectors are 0.8mm appart. Can interpolate detectors 

for y = 1:geo.numDet
 
  detY(count) = detSpace/2 + (geo.detBayHeight/2 - (y*detSpace));
  count = count + 1;
  
end

%Define block face.  radFace = 0 means flat face
switch(geo.lensType) 
  case 'circle'
    hFace = (2*(-d1) + ((2*(-d1))^2 - 4*((-d1)^2 - geo.radFace^2))^0.5)/2;
    kFace = 0;
  case 'ellipse'
    hFace = (-d1)+(geo.bEll^2/(1 + geo.ecc^2))^0.5; 
    kFace = 0;
    
    aEll = (geo.bEll^2/(1 + geo.ecc^2))^0.5;
end

%Define REAR block face.
    hFace2 = (d2)-(geo.bEll2^2/(1 + geo.ecc2^2))^0.5; 
    kFace2 = 0;
    
    aEll2 = (geo.bEll2^2/(1 + geo.ecc2^2))^0.5;


%################# Compute Intersection Matrices ##########################

% define size of intersection matrices
% XintersectionMatrix = nan(r,10);
% YintersectionMatrix = nan(r,10);
XintersectionMatrix = zeros(r,10);
YintersectionMatrix = zeros(r,10);

% loop over rays
for i = 1:r
  
  % define laser source
  switch(geo.laserType)
     case 'optimalCrossover'
      XintersectionMatrix(i,10) = x0;
      YintersectionMatrix(i,10) = y0;
    case 'fan'
      XintersectionMatrix(i,10) = x0;
      YintersectionMatrix(i,10) = y0;
    case 'parallel'
      x0 = geo.x0(i);
      y0 = geo.y0(i);
      
      XintersectionMatrix(i,10) = x0;
      YintersectionMatrix(i,10) = y0;
      
      angles = zeros(1,r);
    case('crossover')
      x0 = geo.x0(i);
      y0 = geo.y0(i);
      
      XintersectionMatrix(i,10) = x0;
      YintersectionMatrix(i,10) = y0;
      
      
      angles = ones(1,r) * atan(((detY(geo.bottomDetYcross)-(detSpace/2)) - geo.y0(1))/(det - geo.x0(1)));
  end
  
  
  %ray describes the parametertes of the ray [x y slope]
  ray = [x0 y0 tan(angles(i))];
  
  switch(geo.lensType) 
    case 'circle'
      if geo.radFace == 0
        %%%%%% Flat face
        % define intersection points with front edge of acrylic block
        XintersectionMatrix(i,9) = wall;
        YintersectionMatrix(i,9) = ray(3)*(wall - ray(1)) + ray(2);

        % REFRACTION (Air -> AcrylicBlock)
        % call the Snell's Law function to calculate the angle of refraction
        aRef = Snells(air,acrylic,angles(i));

        % define parameters of new ray after refraction 
        ray = [wall (ray(3)*(wall - ray(1)) + ray(2)) tan(aRef)];

        % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
        %pick left or right intersection
        xint = xintl;
        yint = yintl;

      else
        %%%%%% Circle face
        % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),hFace,kFace,geo.radFace);
        %pick left or right intersection
        xint = xintl;
        yint = yintl;

        if discrim > 0

          % define intersection points with front edge of acrylic block
          XintersectionMatrix(i,9) = xint;
          YintersectionMatrix(i,9) = yint;

          % calculate angle of incidence wrt normal of Air->Acrylic interface
          aInc = atan((tan(angles(i)) - (yint-kFace)/(xint-hFace))/(1+(yint-kFace)/(xint-hFace)*tan(angles(i))));

          % REFRACTION (Air -> AcrylicBlock)
          % call the Snell's Law function to calculate the angle of refraction
          %           aRef = -(aInc - angles(i) - Snells(air,acrylic,aInc));
          aRef = Snells(air,acrylic,aInc)+atan((yint-kFace)/(xint-hFace));
          ray = [xint yint tan(aRef)];

          % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
          [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
          %pick left or right intersection
          xint = xintl;
          yint = yintl;
        end
      end
      
    case 'ellipse'
      [xintl,yintl,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace,kFace,geo.ecc,geo.bEll);
      %pick left or right intersection
      xint = xintl;
      yint = yintl;

      if discrim > 0

        % define intersection points with front edge of acrylic block
        XintersectionMatrix(i,9) = xint;
        YintersectionMatrix(i,9) = yint;

        % calculate angle of incidence wrt normal of air->acrylic interface
        aInc = atan((tan(angles(i)) - (((aEll^2)*(yint-kFace))/((geo.bEll^2)*(xint-hFace)))) / (1 + (((aEll^2)*(yint-kFace))/((geo.bEll^2)*(xint-hFace)))*tan(angles(i))));

        % REFRACTION (Air -> AcrylicBlock)
        % call the Snell's Law function to calculate the angle of refraction
        aRef = -(aInc - angles(i) - (Snells(air,acrylic,aInc,geo.polAng)));
        ray = [xint yint tan(aRef)];

        % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
        %pick left or right intersection
        xint = xintl;
        yint = yintl;
      end
  end

  % seperate rays that do not intersect Acrylic->Bath interface
  if discrim > 0
    
    % define intersection points with Acrylic->Bath interface
    XintersectionMatrix(i,8) = xint;
    YintersectionMatrix(i,8) = yint;
    
    % calculate angle of incidence wrt normal of Acrylic->Bath interface
    aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
    
    % determain if rays are reflected or refracted
    if abs(aInc) < crit1
      
      % REFRACTION (AcrylicBlock -> Bath)
      aRef = -(aInc - atan(ray(3)) - Snells(acrylic,bath,aInc,geo.polAng));
      ray = [xint yint tan(aRef)];
      % calculate intersection points with Bath->Container interface
      [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
      xint = xintl;
      yint = yintl;
      
      % seperate rays that do not intersect Bath->Container interface
      if discrim > 0
        
         % define intersection points with Bath->Container Interface 
         XintersectionMatrix(i,7) = xint;
         YintersectionMatrix(i,7) = yint;
         % calculate angle of incidence wrt normal of Bath->Container interface
         aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
         
           
         % REFRACTION (Bath -> Container)
         aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc,geo.polAng));
         ray = [xint yint tan(aRef)];
         %calculate intersection points with Container->Gel interface
         [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
         xint = xintl;
         yint = yintl;

         % seperate rays that do not intersect Container->Gel interface
         if discrim > 0 
           
           % define intersection points with Container->Gel Interface 
             XintersectionMatrix(i,6) = xint;
             YintersectionMatrix(i,6) = yint;
             % calculate angle of incidence wrt normal of Container->Gel interface
             aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
           
           if abs(aInc) < crit2
             
             % REFRACTION (Container -> Gel)
             aRef = -(aInc - atan(ray(3)) - Snells(acrylic,gel,aInc,geo.polAng));
             ray = [xint yint tan(aRef)];
             %calculate intersection points with Gel->ContainerBack interface
             [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
             xint = xintr;
             yint = yintr;
             % define intersection points with Gel->ContainerBack
             XintersectionMatrix(i,5) = xint;
             YintersectionMatrix(i,5) = yint;
             % calculate angle of incidence wrt normal of Gel->ContainerBack interface
             aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
             
             % REFRACTION (Gel -> ContainerBack)
             aRef = -(aInc - atan(ray(3)) - Snells(gel,acrylic,aInc,geo.polAng));
             ray = [xint yint tan(aRef)];
             %calculate intersection points with Container->BathBack interface
             [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
             xint = xintr;
             yint = yintr;
             % define intersection points with Container->BathBack
             XintersectionMatrix(i,4) = xint;
             YintersectionMatrix(i,4) = yint;
             % calculate angle of incidence wrt normal of Container->BathBack interface
             aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
               
             if abs(aInc) < crit1
               
               % REFRACTION (Container -> BathBack)
               aRef = -(aInc - atan(ray(3)) - Snells(acrylic,bath,aInc,geo.polAng));
               ray = [xint yint tan(aRef)];
               %calculate intersection points with Bath->AcrylicBack interface
               [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
               xint = xintr;
               yint = yintr;
               % define intersection points with Bath->AcrylicBack
               XintersectionMatrix(i,3) = xint;
               YintersectionMatrix(i,3) = yint;
               % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
               aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
               
               % REFRACTION (Bath -> AcrylicBack)
               aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc,geo.polAng));
               ray = [xint yint tan(aRef)];
%                % define intersection points with Detectors
%                XintersectionMatrix(i,2) = det; %supposed to be rear wall
%                YintersectionMatrix(i,2) = ray(3)*(det - ray(1)) + ray(2);
               
               [xintl,yintl,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace2,kFace2,geo.ecc2,geo.bEll2);
               %pick left or right intersection
               xint = xintr;
               yint = yintr;
               
               if discrim > 0
                   
                   % define intersection points with back edge of acrylic block
                   XintersectionMatrix(i,2) = xint;
                   YintersectionMatrix(i,2) = yint;
                   
                   % calculate angle of incidence wrt normal of acrylic->air interface
                   aInc = atan((tan(angles(i)) - (((aEll2^2)*(yint-kFace2))/((geo.bEll2^2)*(xint-hFace2)))) / (1 + (((aEll2^2)*(yint-kFace2))/((geo.bEll2^2)*(xint-hFace2)))*tan(angles(i))));
                   
                   
                   if abs(aInc) < crit3 %%%%%%%%%new
                       
                   % REFRACTION (block -> air)
                   aRef = -(aInc - atan(ray(3)) - Snells(acrylic,air,aInc,geo.polAng));
                   ray = [xint yint tan(aRef)];
                   %calculate intersection points with Bath->AcrylicBack interface
                   % define intersection points with Detectors
                   XintersectionMatrix(i,1) = det; %supposed to be rear wall
                   YintersectionMatrix(i,1) = ray(3)*(det - ray(1)) + ray(2);
                   else
                       %Totaly internally reflected rays (Acrylic -> Air)
                       
                       % REFLECTION (block -> air)
                       aRef = pi - aInc - (aInc - atan(ray(3)));
                       ray = [xint yint tan(aRef)];
                       if aRef > pi
                           aRef = aRef-2*pi;
                       end
                       
                       if aRef > -pi/2 && aRef < pi/2
                           [xintl,yintl,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace2,kFace2,ecc2,bEll2);
                           %pick left or right intersection
                           xint = xintr;
                           yint = yintr;
                           
                           if discrim > 0
                               
                               % define intersection points with back edge of acrylic block
                               XintersectionMatrix(i,1) = xint;
                               YintersectionMatrix(i,1) = yint;
                           end
                           
                           
                       end
                   end
               end
               
             else
               %Totaly internally reflected rays (Container -> BathBack)
               
               
             end
             
           else
             %Totaly internally reflected rays (Container -> Gel)
             
             % REFLECTION (Container -> Gel)
             aRef = pi - aInc - (aInc - atan(ray(3)));
             ray = [xint yint tan(aRef)];

             %calculate intersection points with Container->BathBack interface
             [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
             xint = xintr;
             yint = yintr;
             % define intersection points with Container->BathBack
             XintersectionMatrix(i,5) = xint;
             YintersectionMatrix(i,5) = yint;
             % calculate angle of incidence wrt normal of Container->BathBack interface
             aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
               
             if abs(aInc) < crit1
               
               % REFRACTION (Container -> BathBack)
               aRef = -(aInc - atan(ray(3)) - Snells(acrylic,bath,aInc,polAng));
               ray = [xint yint tan(aRef)];
               %calculate intersection points with Bath->AcrylicBack interface
               [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
               xint = xintr;
               yint = yintr;
               % define intersection points with Bath->AcrylicBack
               XintersectionMatrix(i,4) = xint;
               YintersectionMatrix(i,4) = yint;
               % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
               aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
               
               % REFRACTION (Bath -> AcrylicBack)
               aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc,polAng));
               ray = [xint yint tan(aRef)];
               % define intersection points with Detectors
               XintersectionMatrix(i,3) = det;
               YintersectionMatrix(i,3) = ray(3)*(det - ray(1)) + ray(2);
               
             else
               %Totaly internally reflected rays (Container -> BathBack)
               % *Ignore* - will not occur with any geometry
             end
             
           end
           
         % rays that go straight through Container
         else
           
           % calculate intersection points with Container->BathBack
           [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
           xint = xintr;
           yint = yintr;
           % define intersection points of rays that go straight through Container
           XintersectionMatrix(i,6) = xint;
           YintersectionMatrix(i,6) = yint;
           % calculate angle of incidence wrt normal of Container->BathBack interface
           aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
           
           if abs(aInc) < crit1 
             
             %REFRACTION (Container->BathBack)
             aRef = -(aInc - atan(ray(3)) - Snells(acrylic,bath,aInc));
             ray = [xint yint tan(aRef)];
             %calculate intersection points with Bath->AcrylicBack interface
             [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
             xint = xintr;
             yint = yintr;
             % define intersection points with Bath->AcrylicBack
             XintersectionMatrix(i,5) = xint;
             YintersectionMatrix(i,5) = yint;
             % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
             aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
             
             %REFRACTION (Bath->AcrylicBack)
             aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc));
             ray = [xint yint tan(aRef)];
             % define intersection points with Detectors
             XintersectionMatrix(i,4) = det;
             YintersectionMatrix(i,4) = ray(3)*(det - ray(1)) + ray(2);
             
           else
             %Totaly internanly reflected rays (Container->BathBack)
             % Rays that go straight through container wall
             % *Ignore* - will not occur with any geometry
           end
              
         end
           
      % rays that go straight through Bath  
      else
        
        % calculate intersection points with Bath->AcrylicBack
        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
        xint = xintr;
        yint = yintr;
        % define intersection points of rays that go straight through Bath
        XintersectionMatrix(i,7) = xint;
        YintersectionMatrix(i,7) = yint;
        % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
        
        % REFRACTION (Bath -> AcrylicBack)
        aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc,geo.polAng));
        ray = [xint yint tan(aRef)];
        % define intersection points with Detectors
        XintersectionMatrix(i,6) = det;
        YintersectionMatrix(i,6) = ray(3)*(det - ray(1)) + ray(2);
        
      end
       
    else
      %Totaly internaly reflected rays (AcrylicBlock -> Bath)
      
      % REFLECTION (AcrylicBlock -> Bath)
      aRef = pi - aInc - (aInc - atan(ray(3)));
      ray = [xint yint tan(aRef)];
     
      %calculate intersection points with Bath->AcrylicBack interface
      [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
      xint = xintr;
      yint = yintr;
      % define intersection points with Bath->AcrylicBack
      XintersectionMatrix(i,7) = xint;
      YintersectionMatrix(i,7) = yint;
      % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
      aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

      % REFRACTION (Bath -> AcrylicBack)
      aRef = -(aInc - atan(ray(3)) - Snells(bath,acrylic,aInc,geo.polAng));
      ray = [xint yint tan(aRef)];
      % define intersection points with Detectors
      XintersectionMatrix(i,6) = det;
      YintersectionMatrix(i,6) = ray(3)*(det - ray(1)) + ray(2);
               
    end
    
  % rays that go straight to Detectors
  else
    
    % define intersection points of rays thats go straight to detectors
    XintersectionMatrix(i,8) = det;
    YintersectionMatrix(i,8) = ray(3)*(det - ray(1)) + ray(2);
    
  end
  
end

%Extract only numDet rays (ray closest to middle of each detector)

% %Create indice of rays that reach the detectors
% I = [];
% J = [];
% for i = 1:size(XintersectionMatrix,1)
%   for j = 1:size(XintersectionMatrix,2)
%     
%     if XintersectionMatrix(i,j) == det && YintersectionMatrix(i,j) <= detY(1)...
%         && YintersectionMatrix(i,j) >= detY(end)
%       
%       I = [I i];
%       J = [J j];
%       
%     end
%     
%   end
% end
% 
% %Select rays closest to the middle of each detctor
% rayNumber = [];
% 
% for l = 1:length(detY)
%   
%   minVal = inf;
%   rayNum = [];
%   for p = 1:length(I)
% 
%     if abs(YintersectionMatrix(I(p),J(p))-detY(l)) < minVal
% 
%          minVal = abs(YintersectionMatrix(I(p),J(p))-detY(l));
%          rayNum = I(p);
%     end
%   end
%   rayNumber = [rayNumber rayNum];
% end
% 
% 
% xInts = [];
% yInts = [];
% for a = 1:length(rayNumber)
%   
%   xInts = [xInts ; XintersectionMatrix(rayNumber(a),:)];
%   yInts = [yInts ; YintersectionMatrix(rayNumber(a),:)];
%   
% end
% 
% dupRayIndx = [];
% for i = 1:size(xInts,1)-1
%   
%   if isequaln(xInts(i,:),xInts(i+1,:)) && isequaln(yInts(i,:),yInts(i+1,:))
% 
%     dupRayIndx = [dupRayIndx i];
%     
%   end
%   
% end
%   
% xInts(dupRayIndx,:) = nan;
% yInts(dupRayIndx,:) = nan;
% 
% assignin('base', 'xInts', xInts)
% assignin('base', 'yInts', yInts)
%   
% tEnd = toc(tStart);
% fprintf('Complete Refraction Function Is Detected Zero Duplicates Time -> %f sec\n', tEnd);
% 


return
end




