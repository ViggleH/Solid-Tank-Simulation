%##########################################################################
%                     Image Making Function
%##########################################################################

function [I] = Image(type,size,IntMatrix,i,ray,imCoordx,imCoordy,scale,imangles,xint,yint,r1,theta,h,k,u);
%Discription
%   This function takes inputs for the type of image, along with size and
%   location of the image. it then determines if this ray crosses one or
%   more of the images, find the distance across and then attenuates the intesity 
%   appropatly
%attenuation coefficent in inverse millimeters
%u = 0.0166; %Flexidos3d
% u = 0.0110; %ClearView
% u = 0.0; %Water

%% circle Type 1

if type == 1
    
    [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),imCoordx,imCoordy,size);%size gives the radius of the circle
    
    if discrim>0 %determines if rays intersect circle then changes the intesity matrix
        distance=sqrt((xintr-xintl)^2+(yintr-yintl)^2)/scale;%calculates the distance the ray crosses inside the circle
        I=IntMatrix*exp(-u*distance);
        
        
        
        %Plotting the circle
        %figure(1)
        %hold on
        %plot(imCoordx+r4*sin(theta),imCoordy+r4*cos(theta), 'LineWidth',1, 'Color', 'red');
    else %intesity remains the same
        I = IntMatrix;
    end
end

%% ellipse Type 2

if type == 2;
  
%defining top and bottom "points/corners" of ellipse
    length=20.692;
    radiusToImage=20*scale;
    distance=sqrt((radiusToImage+(sqrt(3)*length/2))^2+(length/2)^2);
    theta=sin((length/2)/(distance));
    theta2=imangles+theta;
    theta3=imangles-theta;
%"top" point
    x1=distance*cos(theta2);
    y1=distance*sin(theta2);
%"bottom" point
    x2=distance*cos(theta3);
    y2=distance*sin(theta3);
    
%Data for ellipse
    e=0.99; %eccentricity of the elipse
    dia = 104; %diameter of bore hole
    d1 = 183.5+(dia/2); % distance from front edge of acrylic block to center of bore
    a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
    b = a*sqrt(1-e^2);
    t = 0:(pi/2000):2*pi;
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(y2-y1,x2-x1);
    xbox = ((x1+x2)/2 + X*cos(w) - Y*sin(w));%Defines the elipse
    ybox = ((y1+y2)/2 + X*sin(w) + Y*cos(w));%Defines the elipse
    
figure(1)
hold on
plot(xbox,ybox,'Color', 'red');

%finding intersecting lines  
    [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
    xx = [xint xintr];
    yy = [yint yintr];
    [xi,yi] = polyxpoly(xx,yy,xbox,ybox);%this part takes lots of the time 
    

    xi=[xi xi];
    yi=[yi yi];
    if xi~=0%determine is the ray passes through the ellipse then find attenuation
        distance=sqrt((xi(1,1)-xi(2,2))^2+(yi(1,1)-yi(2,2))^2)/scale;%find the distance between the point of entry and point of leaving 
        I=IntMatrix*exp(-u*distance);
        IntMatrix=I;
    end

% figure(1)
% hold on
% plot(xbox,ybox,'Color', 'red');
%mapshow(xi,yi,'DisplayType','point','Marker','o')
 end

%% square Type 3


if type == 3 ;
    
%defining all four corners of square
    length=15.692;%defines length of one side of square
    radiusToImage=20;%distance from center of gel to closest point of square
    distance=sqrt((radiusToImage+(sqrt(3)*length/2))^2+(length/2)^2);
    theta=sin((length/2)/(distance));
    theta2=imangles+theta;
    theta3=imangles-theta;
%closest two points
    X1=distance*cos(theta2);
    Y1=distance*sin(theta2);
    X2=distance*cos(theta3);
    Y2=distance*sin(theta3);
%furthest two points
    X3=X2+cos(imangles)*length;
    Y3=Y2+sin(imangles)*length;
    X4=X1+cos(imangles)*length;
    Y4=Y1+sin(imangles)*length;
%putting points in matrix defining the square        
    xbox=[X1,X2,X3,X4,X1];
    ybox=[Y1,Y2,Y3,Y4,Y1];
%finding intersecting lines                                    
    [xintl,yintl,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
    xx = [xint xintr];
    yy = [yint yintr];
    [xi,yi] = polyxpoly(xx,yy,xbox,ybox);%this part takes lots of the time 
                                  
    
%     figure(1)
%     hold on
%     plot(xbox,ybox,'Color', 'red')
%     mapshow(xi,yi,'DisplayType','point','Marker','o')

    xi=[xi xi];
    yi=[yi yi];
    
    if xi~=0%determine is the ray passes through the square then find attenuation
        distance=sqrt((xi(1,1)-xi(2,2))^2+(yi(1,1)-yi(2,2))^2)/scale;
        I=IntMatrix*exp(-u*distance);
        IntMatrix=I;
    end
    
    
end