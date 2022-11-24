%##########################################################################
%                     COMPLETE REFRACTION FUNCTION
%##########################################################################

function [XintersectionMatrix,YintersectionMatrix,IntMatrix,angMatrix] = ...
    RefracFunc_Merge(N,r,gel,span1,span2,scale,dlaser,bath,dia,d1,d2, ...
    lensType,ecc,bEll,h,k,radFace,polAng,imCoordx,imCoordy,imangles,theta,r4, u)

%######################### Function Inputs ################################

% N = N; % image size
% r = r; % number of rays
% gel = gel; % refractive index of the gel
% span1 = span1; % the span of the fan beam(in radian)
% span2 = span2; % the span of the fan beam(in radian)
% scale = scale; % scale factor relating the physical paramters of the syetem to the image size
% dlaser = dlaser; % distance from laser source to acrylic block
% bath = bath; % refractive index of the matching bath
% radFace = 
% bEll = 

%################## Physical parameters of system #########################

% dia = 104; % 104mm - diameter of bore hole
wallT = 0.125*25.4; %0.125in wall - container wall thickness
gapT = (103.6-dia)/2; % gap thickness i.e. bath thickness
% d1 = 183.5+(dia/2); % distance from front edge of acrylic block to center of bore
% d2 = 344-d1; % distance from center of bore hole to detectors - 344mm = total length of acrylic block
doff = 0; % vertical offset of laser source
r0 = 1; %defines the initial ray intesity to be one

h = h*scale; k = k*scale; % center of bore at origin
x0 = -d1*scale - dlaser*scale; % laser source position on x-axis
det = d2*scale; % detectors position
y0 = -doff*scale; % laser source position on y-axis
wall = -d1*scale; % position of front edge of acrylic block

angles = linspace(span1,span2,r); % angle of inclination of each ray from laser source
r1 = (dia/2)*scale - wallT*scale; %gel radius
r2 = (dia/2)*scale;  %container radius
r3 = (dia/2)*scale + gapT*scale; %bore radius


air = 1.00029; % refractive index of air
acrylic = 1.4886; % refractive index of acrylic (block and container)

% critical angles assuming bath < acrylic
crit1 = asin(bath/acrylic); % critical angle where total internal reflection occurs at Acrylic->Bath interface
crit2 = asin(gel/acrylic); % critical angle where total internal reflection occurs at Gel->Acrylic interface

aEll = (bEll^2/(1 + ecc^2))^0.5;

%define curved face.  radius = 0 means flat face
x = -d1*scale;
if lensType == 1 %circle
    hFace = (2*x + ((2*x)^2 - 4*(x^2 - radFace^2))^0.5)/2;
    kFace = 0;
elseif lensType == 2 %ellipse
    hFace = x+aEll; %b term from ellipse equation using a and ecc
    kFace = 0;
end

%################# Compute Intersection Matrices ##########################

% define size of intersection matrices
XintersectionMatrix = zeros(r,9);
YintersectionMatrix = zeros(r,9);
IntMatrix = zeros(r,9);
angMatrix = zeros(r,9);

% loop over rays
for i = 1:r
    
    % define laser source
    XintersectionMatrix(i,9) = x0;
    YintersectionMatrix(i,9) = y0;
    IntMatrix(i,9) = real(r0);
    angMatrix(i,9) = real(angles(i)*(180/pi));
    
    %ray describes the parametertes of the ray [x y slope]
    ray = [x0 y0 tan(angles(i))];
    if lensType == 1
        if radFace == 0
            %%%%%% Flat face
            % define intersection points with front edge of acrylic block
            XintersectionMatrix(i,8) = wall;
            YintersectionMatrix(i,8) = ray(3)*(wall - ray(1)) + ray(2);
            
            % REFRACTION (Air -> AcrylicBlock)
            % call the Snell's Law function to calculate the angle of refraction
            [aRef,iRef] = Snells(air,acrylic,angles(i),polAng);
            angMatrix(i,8)=aRef*(180/pi);
            IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
            
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
            [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),hFace,kFace,radFace);
            %pick left or right intersection
            xint = xintl;
            yint = yintl;    % seperate rays that do not intersect with curved block face
            
            if discrim > 0
                
                % define intersection points with front edge of acrylic block
                XintersectionMatrix(i,8) = xint;
                YintersectionMatrix(i,8) = yint;
                
                % calculate angle of incidence wrt normal of air->acrylic interface
                aInc = atan((tan(angles(i)) - (yint-kFace)/(xint-hFace))/(1+(yint-kFace)/(xint-hFace)*tan(angles(i))));
                
                % same expression as above but smaller.This is the Andy version
                % aInc = atan(ray(3)) + atan((yint - kFace)/(xint - hFace))
                
                % REFRACTION (Air -> AcrylicBlock)
                % call the Snell's Law function to calculate the angle of refraction
                % aRef = Snells(air,acrylic,aInc);
                % ray = [xint yint tan(aRef--atan((yint-kFace)/(xint-hFace)))];
                [p,iRef] = Snells(air,acrylic,aInc,polAng);
                aRef = p+atan((yint-kFace)/(xint-hFace));
                angMatrix(i,8)=aRef*(180/pi);
                IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
                ray = [xint yint tan(aRef)];
                
                % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                %pick left or right intersection
                xint = xintl;
                yint = yintl;
            end
        end
    elseif lensType == 2
        [xintl,yintl,xintr,yintr,discrim] = Lens(ray(1),ray(2),ray(3),hFace,kFace,ecc,bEll);
        %pick left or right intersection
        xint = xintl;
        yint = yintl;
        
        if discrim > 0
            
            % define intersection points with front edge of acrylic block
            XintersectionMatrix(i,8) = xint;
            YintersectionMatrix(i,8) = yint;
            
            % calculate angle of incidence wrt normal of air->acrylic interface
%             aInc = atan((tan(angles(i)) - (yint-kFace)/(xint-hFace))/(1+(yint-kFace)/(xint-hFace)*tan(angles(i))));
            aInc = atan((tan(angles(i)) - (((aEll^2)*(yint-kFace))/((bEll^2)*(xint-hFace)))) / (1 + (((aEll^2)*(yint-kFace))/((bEll^2)*(xint-hFace)))*tan(angles(i))));
            % REFRACTION (Air -> AcrylicBlock)
            % call the Snell's Law function to calculate the angle of refraction
            % aRef = Snells(air,acrylic,aInc);
            % ray = [xint yint tan(aRef--atan((yint-kFace)/(xint-hFace)))];
            [p,iRef] = Snells(air,acrylic,aInc,polAng);
            aRef = p+atan((yint-kFace)/(xint-hFace));
            angMatrix(i,8)=aRef*(180/pi);
            IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
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
        XintersectionMatrix(i,7) = xint;
        YintersectionMatrix(i,7) = yint;
                
        % calculate angle of incidence wrt normal of Acrylic->Bath interface
        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
        
        % determain if rays are reflected or refracted
        if abs(aInc) < crit1
            
            % REFRACTION (AcrylicBlock -> Bath)
            [p,iRef] = Snells(acrylic,bath,aInc,polAng);
            aRef = -(aInc - aRef - p);
            angMatrix(i,7)=aRef*(180/pi);
            IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
            ray = [xint yint tan(aRef)];
            % calculate intersection points with Bath->Container interface
            [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
            xint = xintl;
            yint = yintl;
            
            % seperate rays that do not intersect Bath->Container interface
            if discrim > 0
                
                % define intersection points with Bath->Container Interface
                XintersectionMatrix(i,6) = xint;
                YintersectionMatrix(i,6) = yint;
                
                % calculate angle of incidence wrt normal of Bath->Container interface
                aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                
                
                % REFRACTION (Bath -> Container)
                [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                aRef = -(aInc - aRef - p);
                angMatrix(i,6)=aRef*(180/pi);
                IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
                ray = [xint yint tan(aRef)];
                %calculate intersection points with Container->Gel interface
                [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
                xint = xintl;
                yint = yintl;
                
                % seperate rays that do not intersect Container->Gel interface
                if discrim > 0
                    
                    % define intersection points with Container->Gel Interface
                    XintersectionMatrix(i,5) = xint;
                    YintersectionMatrix(i,5) = yint;
                    
                    % calculate angle of incidence wrt normal of Container->Gel interface
                    aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                    
                    if abs(aInc) < crit2
                        
                        % REFRACTION (Container -> Gel)
                        [p,iRef] = Snells(acrylic,gel,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,5)=aRef*(180/pi);
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        
%%%%%%%%%%%%% new part finding lines crossing circle and degrading there intesity accourdingly                         
                        %calculate intersection points with image  
                        %type 1=circle
                        %type 2 = ellipse
                        %type 3 = square
%                         type=3;
%                         [IntM] = Image(type,r4,IntMatrix(i,5),i,ray,imCoordx,imCoordy,scale,imangles,xint,yint,aRef,det,r1,theta);
%                         IntMatrix(i,5)=IntM;
%                         
%                          type=2;
%                         [IntM] = Image(type,r4,IntMatrix(i,5),i,ray,imCoordx,imCoordy,scale,imangles,xint,yint,aRef,det,r1,theta);
%                         IntMatrix(i,5)=IntM;
%                         

%flexidose image
                        type=1;
                        [IntM] = Image(type,r4,IntMatrix(i,5),i,ray,imCoordx,imCoordy,scale,imangles,xint,yint,aRef,det,r1,theta,h,k, u);
                        IntMatrix(i,5)=IntM;
                        
                        %calculate intersection points with Gel->ContainerBack interface
                        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Gel->ContainerBack
                        XintersectionMatrix(i,4) = xint;
                        YintersectionMatrix(i,4) = yint;
                        
                        % calculate angle of incidence wrt normal of Gel->ContainerBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                        
                        % REFRACTION (Gel -> ContainerBack)
                        [p,iRef] = Snells(gel,acrylic,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,4)=aRef*(180/pi);
                        IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        %calculate intersection points with Container->BathBack interface
                        [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Container->BathBack
                        XintersectionMatrix(i,3) = xint;
                        YintersectionMatrix(i,3) = yint;
                        
                        % calculate angle of incidence wrt normal of Container->BathBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                        
                        if abs(aInc) < crit1
                            
                            % REFRACTION (Container -> BathBack)
                            [p,iRef] = Snells(acrylic,bath,aInc,polAng);
                            aRef = -(aInc - aRef - p);
                            angMatrix(i,3)=aRef*(180/pi);
                            IntMatrix(i,3)=IntMatrix(i,4)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            %calculate intersection points with Bath->AcrylicBack interface
                            [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                            xint = xintr;
                            yint = yintr;
                            % define intersection points with Bath->AcrylicBack
                            XintersectionMatrix(i,2) = xint;
                            YintersectionMatrix(i,2) = yint;
                            
                            % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                            aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                            
                            % REFRACTION (Bath -> AcrylicBack)
                            [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                            aRef = -(aInc - aRef - p);
                            angMatrix(i,2)=aRef*(180/pi);
                            IntMatrix(i,2)=IntMatrix(i,3)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            % define intersection points with Detectors
                            XintersectionMatrix(i,1) = det;
                            YintersectionMatrix(i,1) = ray(3)*(det - ray(1)) + ray(2);
%                             angMatrix(i,1)=aRef*(180/pi);
%                             IntMatrix(i,1)=IntMatrix(i,2)*(1-iRef);
                        else
                        end
                        
                    else
                    end
                    
                    % rays that go straight through Container
                else
                    
                    % calculate intersection points with Container->BathBack
                    [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
                    xint = xintr;
                    yint = yintr;
                    % define intersection points of rays that go straight through Container
                    XintersectionMatrix(i,5) = xint;
                    YintersectionMatrix(i,5) = yint;
                    
                    % calculate angle of incidence wrt normal of Container->BathBack interface
                    aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                    
                    if abs(aInc) < crit1
                        
                        %REFRACTION (Container->BathBack)
                        [p,iRef] = Snells(acrylic,bath,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,5)=aRef*(180/pi);
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
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
                        
                        %REFRACTION (Bath->AcrylicBack)
                        [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                        aRef = -(aInc - aRef - p);
                        angMatrix(i,4)=aRef*(180/pi);
                        IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        % define intersection points with Detectors
                        XintersectionMatrix(i,3) = det;
                        YintersectionMatrix(i,3) = ray(3)*(det - ray(1)) + ray(2);
%                         angMatrix(i,3)=aRef*(180/pi);
%                         IntMatrix(i,3)=IntMatrix(i,4)*(1-iRef);
                        
                    else
                    end
                    
                end
                
                % rays that go straight through Bath
            else
                
                % calculate intersection points with Bath->AcrylicBack
                [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                xint = xintr;
                yint = yintr;
                % define intersection points of rays that go straight through Bath
                XintersectionMatrix(i,6) = xint;
                YintersectionMatrix(i,6) = yint;
                
                % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));
                
                % REFRACTION (Bath -> AcrylicBack)
                [p,iRef] = Snells(bath,acrylic,aInc,polAng);
                aRef = -(aInc - aRef - p);
                angMatrix(i,6)=aRef*(180/pi);
                IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
                ray = [xint yint tan(aRef)];
                % define intersection points with Detectors
                XintersectionMatrix(i,5) = det;
                YintersectionMatrix(i,5) = ray(3)*(det - ray(1)) + ray(2);
%                 angMatrix(i,5)=aRef*(180/pi);
%                 IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
                
            end
            
        else
        end
        
        % rays that go straight to Detectors
    else
        
        % define intersection points of rays thats go straight to detectors
        XintersectionMatrix(i,7) = det;
        YintersectionMatrix(i,7) = ray(3)*(det - ray(1)) + ray(2);
%         angMatrix(i,7)=aRef*(180/pi);
%         IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
    end
end
%         end
%     end
% end
return