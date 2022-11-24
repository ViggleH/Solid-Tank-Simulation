%##########################################################################
%          COMPLETE REFRACTION FUNCTION IS DETECTED ZERO DUPLICATES
%##########################################################################

function [XintersectionMatrix,YintersectionMatrix,intensityMat,angleMat] = RefracFunc_CCD_Sim(x, gel, u)

tStart = tic;
%################## Physical parameters of system #########################
Solid_Tank_Sim_Constants

L = x(1);
h = x(2);
dlaser = x(3);
bEll = x(4);
ecc = x(5);
bEll2 = x(6);
ecc2 = x(7);
d3 = x(8);

d1 = L/2 + h;
d2 = L/2 - h;

h = 0;

imCoordx= radiusToImage*cos(imangles)+h*scale;%starting x coord
imCoordy= radiusToImage*sin(imangles)+k*scale;%starting y coord

det = d2 + d3; %need to figure out how to assign this

r0 = 1; %defines the initial ray intesity to be one

switch(laserType)
    case('optimalCrossover')
        x0 = -d1 - dlaser; % laser source position on x-axis
        y0 = doff; % laser source position on y-axis
    case('fan')
        x0 = -d1 - dlaser; % laser source position on x-axis
        y0 = doff; % laser source position on y-axis
    case('parallel')
        x0 = ones(1,r) * (-d1 - dlaser);
        y0 = linspace((doff) - (dia/2) , (doff) + (dia/2) ,r);
    case('line')
        x0 = -d1 - dlaser; % laser source position on x-axis
        y0 = doff; % laser source position on y-axis
end

wall = -d1; % position of front edge of acrylic block

switch(laserType)
    case 'fan'
        angles = linspace(span1,span2,r); % angle of inclination of each ray from laser source
    case 'optimalCrossover'
        angles = [linspace(-40*pi/180,span2,r/2) linspace(span1,40*pi/180,r/2)]; % angle of inclination of each ray from laser source
    case 'line' %line as opposed to fan
        opp = linspace((d1+d2)*tan(span1),(d1+d2)*tan(span2),r);
        adj = linspace((d1+d2),(d1+d2),length(opp));
        denom = (opp.^2 + adj.^2).^0.5;
        angles = asin(opp./denom);
end

% critical angles assuming bath < acrylic
crit1 = asin(bath/acrylic); % critical angle where total internal reflection occurs at Acrylic->Bath interface
crit2 = asin(gel/acrylic); % critical angle where total internal reflection occurs at Acrylic->Gel interface
crit3 = asin(air/acrylic); % critical angle where total internal reflection occurs at Acrylic->air interface

%Define block face.  radFace = 0 means flat face
switch(lensType)
    case 'circle'
        hFace = (2*(-d1) + ((2*(-d1))^2 - 4*((-d1)^2 - radFace^2))^0.5)/2;
        kFace = 0;
    case 'ellipse'
        hFace = (-d1)+(bEll^2/(1 + ecc^2))^0.5;
        kFace = 0;
        
        aEll = (bEll^2/(1 + ecc^2))^0.5;
end

%Define REAR block face.
hFace2 = (d2)-(bEll2^2/(1 + ecc2^2))^0.5;
kFace2 = 0;

aEll2 = (bEll2^2/(1 + ecc2^2))^0.5;


%################# Compute Intersection Matrices ##########################

% define size of intersection matrices
XintersectionMatrix = nan(r,10);
YintersectionMatrix = nan(r,10);
% XintersectionMatrix = zeros(r,10);
% YintersectionMatrix = zeros(r,10);
IntMatrix = nan(r,10);
angMatrix = nan(r,10);

% loop over rays
for i = 1:r
% for i = 6230:6300
    
    % define laser source
    switch(laserType)
        case 'optimalCrossover'
            XintersectionMatrix(i,10) = x0;
            YintersectionMatrix(i,10) = y0;
        case 'fan'
            XintersectionMatrix(i,10) = x0;
            YintersectionMatrix(i,10) = y0;
        case 'parallel'
            x0 = x0(i);
            y0 = y0(i);
            
            XintersectionMatrix(i,10) = x0;
            YintersectionMatrix(i,10) = y0;
            
            angles = zeros(1,r);
        case 'line'
            XintersectionMatrix(i,10) = x0;
            YintersectionMatrix(i,10) = y0;
    end
    
    IntMatrix(i,10) = real(r0);
    angMatrix(i,10) = real(angles(i)*(180/pi));
    
    %ray describes the parametertes of the ray [x y slope]
    ray = [x0 y0 tan(angles(i))];
    
    switch(lensType)
        case 'circle'
            if geo.radFace == 0
                %%%%%% Flat face
                % define intersection points with front edge of acrylic block
                XintersectionMatrix(i,9) = wall;
                YintersectionMatrix(i,9) = ray(3)*(wall - ray(1)) + ray(2);
                
                % REFRACTION (Air -> AcrylicBlock)
                % call the Snell's Law function to calculate the angle of refraction
                [aRef,iRef] = Snells(air,acrylic,angles(i),polAng);
                angMatrix(i,9)=aRef*(180/pi);
                IntMatrix(i,9)=IntMatrix(i,9)*(1-iRef);
                
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
                    [p,iRef] = Snells(air,acrylic,aInc,polAng);
                    aRef = p+atan((yint-kFace)/(xint-hFace));
                    angMatrix(i,9)=aRef*(180/pi);
                    IntMatrix(i,9)=IntMatrix(i,10)*(1-iRef);
                    ray = [xint yint tan(aRef)];
                    
                    % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                    [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                    %pick left or right intersection
                    xint = xintl;
                    yint = yintl;
                end
            end
            
        case 'ellipse'
            [xintl,yintl,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace,kFace,ecc,bEll);
            %pick left or right intersection
            xint = xintl;
            yint = yintl;
            
            if discrim > 0
                
                % define intersection points with front edge of acrylic block
                XintersectionMatrix(i,9) = xint;
                YintersectionMatrix(i,9) = yint;
                
                % calculate angle of incidence wrt normal of air->acrylic interface
                aInc = atan((tan(angles(i)) - (((aEll^2)*(yint-kFace))/((bEll^2)*(xint-hFace)))) / (1 + (((aEll^2)*(yint-kFace))/((bEll^2)*(xint-hFace)))*tan(angles(i))));
                
                % REFRACTION (Air -> AcrylicBlock)
                % call the Snell's Law function to calculate the angle of refraction
                [p,iRef] = Snells(air,acrylic,aInc,polAng);
                aRef = p+atan((yint-kFace)/(xint-hFace));
                angMatrix(i,9)=aRef*(180/pi);
                IntMatrix(i,9)=IntMatrix(i,10)*(1-iRef);
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
            [p,iRef] = Snells(acrylic,bath,aInc,polAng);
            aRef = -(aInc - aRef - p);
            angMatrix(i,8)=aRef*(180/pi);
            IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
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
                [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                aRef = -(aInc - aRef - p);
                angMatrix(i,7)=aRef*(180/pi);
                IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
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
                        [p,iRef] = Snells(acrylic,gel,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,6)=aRef*(180/pi);
                        IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        
                        %flexidose image
                        type=1;
                        [IntM] = Image(type,r4,IntMatrix(i,6),i,ray,imCoordx,imCoordy,scale,imangles,xint,yint,r1,theta,h,k,u);
                        IntMatrix(i,6)=IntM;
                        
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
                        [p,iRef] = Snells(gel,acrylic,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,5)=aRef*(180/pi);
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
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
                            [p,iRef] = Snells(acrylic,bath,aInc,polAng);
                            aRef = -(aInc - aRef - p);
                            angMatrix(i,4)=aRef*(180/pi);
                            IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
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
                            [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                            aRef = -(aInc - aRef - p);
                            angMatrix(i,3)=aRef*(180/pi);
                            IntMatrix(i,3)=IntMatrix(i,4)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            %                % define intersection points with Detectors
                            %                XintersectionMatrix(i,2) = det; %supposed to be rear wall
                            %                YintersectionMatrix(i,2) = ray(3)*(det - ray(1)) + ray(2);
                            
                            [xintl,yintl,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace2,kFace2,ecc2,bEll2);
                            %pick left or right intersection
                            xint = xintr;
                            yint = yintr;
                            
                            if discrim > 0
                                
                                % define intersection points with back edge of acrylic block
                                XintersectionMatrix(i,2) = xint;
                                YintersectionMatrix(i,2) = yint;
                                
                                % calculate angle of incidence wrt normal of acrylic->air interface
                                aInc = atan((tan(angles(i)) - (((aEll2^2)*(yint-kFace2))/((bEll2^2)*(xint-hFace2)))) / (1 + (((aEll2^2)*(yint-kFace2))/((bEll2^2)*(xint-hFace2)))*tan(angles(i))));
                                
                                
                                if abs(aInc) < crit3 %%%%%%%%%new
                                    
                                    % REFRACTION (block -> air)
                                    [p,iRef] = Snells(acrylic,air,aInc,polAng);
                                    aRef = -(aInc - aRef - p);
                                    angMatrix(i,2)=aRef*(180/pi);
                                    IntMatrix(i,2)=IntMatrix(i,3)*(1-iRef);
                                    ray = [xint yint tan(aRef)];
                                    %calculate intersection points with Bath->AcrylicBack interface
                                    % define intersection points with Detectors
                                    XintersectionMatrix(i,1) = det; %supposed to be rear wall
                                    YintersectionMatrix(i,1) = ray(3)*(det - ray(1)) + ray(2);
                                    angMatrix(i,1)=aRef*(180/pi);
                                    IntMatrix(i,1)=IntMatrix(i,2);
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
                        [p,iRef] = Snells(acrylic,bath,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,6)=aRef*(180/pi);
                        IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
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
                        [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,5)=aRef*(180/pi);
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
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
                [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                aRef = -(aInc - aRef - p);
                angMatrix(i,7)=aRef*(180/pi);
                IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
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
            [p,iRef] = Snells(bath,acrylic,aInc,polAng);
            aRef = -(aInc - aRef - p);
            angMatrix(i,7)=aRef*(180/pi);
            IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
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
IntMatrix(isnan(IntMatrix))=0;
intensityMat = IntMatrix;

angMatrix(isnan(angMatrix))=0;
angleMat = angMatrix;

XintersectionMatrix(isnan(XintersectionMatrix))=0;

YintersectionMatrix(isnan(YintersectionMatrix))=0;

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
tEnd = toc(tStart);
%fprintf('Refraction Function Time -> %f sec\n', tEnd);
%


return
end




