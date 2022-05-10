function gz=gprism(x1,y1,z1,x2,y2,z2,xgv,ygv,zp,rho)
% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% *** from the paper
%%"Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68"
% ****************************************************************

% Analytical solution 
% to compute the vertical attraction of a rectangular prism. 
% Sides of prism are parallel to x,y,z axes, and z is vertical down.
%
% Input parameters:
% Prism extends from x1 to x2, y1 to y2, and z1 to z2
% in x, y, and z directions, respectively.
% xgv: a vector of x coordinates of observation points
% ygv: a vector of y coordinates of observation points
% zp: observation level.
% rho: density of the prism in kg/(m^3).
% All distance parameters in units of km;
%
% Output parameters:
% gz: vertical attraction of the gravity vector g, in mGal.
%
% Reference:
% The algorithm is based on the Subroutine B.6. gbox in the book
% Blakely, R. J., 1996, Potential theory in gravity and magnetic applications:
% Cambridge University Press.
%
gamma=6.67259e-11;
si2mg=1.e5;
km2m=1.e3;

nx=length(xgv);
ny=length(ygv);
gz=zeros(ny,nx);
for ix=1:1:nx
    xp=xgv(ix);
    for jy=1:1:ny
        yp=ygv(jy);
        isign=[-1,1];
        x(1)=xp-x1;
        y(1)=yp-y1;
        z(1)=zp-z1;
        x(2)=xp-x2;
        y(2)=yp-y2;
        z(2)=zp-z2;
        sum=0;
        for i=1:2
            for j=1:2
                for k=1:2
                    rijk=sqrt(x(i)^2+y(j)^2+z(k)^2);
                    ijk=isign(i)*isign(j)*isign(k);
                    arg1=atan((x(i)*y(j))/(z(k)*rijk));
                    if(isnan(arg1))
                        arg1=0;
                    else
                        arg1=z(k)*arg1;
                    end
                    arg2=rijk+y(j);
                    if(arg2==0)
                        arg2=0;
                    else
                        arg2=-x(i)*log(arg2);
                    end
                    arg3=rijk+x(i);
                    if(arg3==0)
                        arg3=0;
                    else
                        arg3=-y(j)*log(arg3);
                    end
                    sum=sum+ijk*(arg1+arg2+arg3);
                end
            end
        end
        gz(jy,ix)=rho*gamma*sum*si2mg*km2m;
    end
end


