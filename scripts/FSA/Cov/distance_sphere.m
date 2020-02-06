function h = distance_sphere(x,y)
%DISTANCE_ Determine distances between locations using the Spherical law of cosines formula
%
%   This function produces a matrix that describes the
%   distances between to sets of locations.
%
%INPUT PARAMETERS
%   x - location coordinates (degrees) for data set #1 [n1 x D] -(long,lat)
%   y - location coordinates (degrees) for data set #2 [n2 x D] -(long,lat)
%OUTPUT PARAMETERS
%   h - distance (km) between points in x from points in y [n1 x n2]
%
%EXAMPLE:
%   x=[0,0; 5,0; 5,5]
%   y=[1,1; 5,5]
%   h = distance_(x,y)
%RESULT:
%    1.4142    7.0711
%    4.1231    5.0000
%    5.6569         0

%EXTERNAL FUNCTIONS CALLED: none
%REVISION HISTORY:
%   pkk, 5/13/97
%   tae, 6/26/98 (minor changes)
%   amm, 5/30/03 this subroutine now uses longitude latitude and calculates distances
%                on a sphere

[n1,D] = size(x);
[n2,D2] = size(y);

if D~=D2
   error('ERROR in DISTANCE_: locations must have same number of dimensions (columns)')
   return;
end


h = zeros(n1,n2);
if D==1
    for id = 1:D
        h = h + (x(:,id)*ones(1, n2)-ones(n1,1)*y(:,id)').^2;
    end
    h = sqrt(h);
else
    r=6371.0087714; %WGS84 mean radius
    x=x*pi/180;
    y=y*pi/180;
    
    for j=1: n2
        
        %h(:,j)=r*acos(sin(x(:,2)).*sin(repmat(y(j,2), n1,1))+cos(x(:,2)).*cos(repmat(y(j,2),n1,1)).*cos(x(:,1)-repmat(y(j,1),n1,1)));
        
        %Transformed GCD to ensure PD
        h(:,j)=2*r*sin(acos(sin(x(:,2)).*sin(repmat(y(j,2), n1,1))+cos(x(:,2)).*cos(repmat(y(j,2),n1,1)).*cos(x(:,1)-repmat(y(j,1),n1,1)))./2);
        
    end
    
    h=real(h);
end
return;

