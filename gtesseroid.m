% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for finding gravity anomalies due to any tesseroid 

function gz=gtesseroid(rad_1,phi_1,lmbda_1,rad_2,phi_2,lmbda_2,Rad,Phi_all,Lmbda_all,rho,y,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %input: 
    %   rad_1       = depth of shallower surface in m
    %   phi_1       = lower limit of longitude of the tesseroid in degree
    %   lmbda_1     = lower limit of latitude of the tesseroid in degree
    %   rad_2       = depth of deeper surface in m
    %   phi_2       = upper limit of longitude of the tesseroid in degree
    %   lmbda_2     = upper limit of latitude of the tesseroid in degree
    %   Rad         = depth of the observation point 
    %   Phi_all     = longitudinal grids of observation points in degree
    %   Lmbda_all   = latitudinal grids of observation points in degree
    %   rho         = density function in kg/m^3
    %   y           = Gaussian quadrature nodes for [-1 1]
    %   w           = Gaussian quadrature weights for [-1 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %output:
    %   gz     = gravity anomaly for given tesseroids at given observations
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %prerequired functions
    %   1. lgwt.m for evaluating nodes and weights 
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Converting all angles from degree to radian
    phi_1=deg2rad(phi_1); phi_2=deg2rad(phi_2);
    lmbda_1=deg2rad(lmbda_1); lmbda_2=deg2rad(lmbda_2);
    Phi_all=deg2rad(Phi_all);  Lmbda_all=deg2rad(Lmbda_all);
    
    %Radius of earth
    R=6371e3;
    %Observation point
    Rad=R-Rad;

    %upper and lower limit of r
    ar=R-rad_1;
    br=R-rad_2;

    %upper and lower limit of phi
    aphi=phi_1;
    bphi=phi_2;

    %upper and lower limit of lambda
    almbd=lmbda_1;
    blmbd=lmbda_2;

    rm=((br-ar)/2).*y+(br+ar)/2;
    
    nx=length(Phi_all);
    ny=length(Lmbda_all);
    gz=zeros(ny,nx);
    %gravitational constant
    G=6.67408*10^-11; %in m^3kg^-1s^-2
    N=length(w);
    for ix=1:1:nx
        Phi=Phi_all(ix);     
        for jy=1:1:ny
            Lmbda=Lmbda_all(jy);
            ss=0;
            %function for which integration have to do
            cosxi=@(phi1,lmbd1) (sin(Phi).*sin(phi1)+cos(Phi).*cos(phi1).*cos(lmbd1-Lmbda));
            f=@(r1,phi1,lmbd1) ((r1.*cosxi(phi1,lmbd1)-Rad)./(sqrt(r1.^2+Rad.^2-2.*Rad.*r1.*cosxi(phi1,lmbd1))).^3).*(r1.^2.*cos(phi1));
            
            for i=1:N
                for j=1:N
                    for k=1:N

                        ym(j)=((bphi-aphi)./2).*y(j)+(bphi+aphi)./2;
                        zm(k)=((blmbd-almbd)./2).*y(k)+(blmbd+almbd)./2;

                        C(i,j,k)=w(i)*w(j)*w(k)*((br-ar)/2).*((bphi-aphi)/2).*((blmbd-almbd)/2);
                        ss=ss+C(i,j,k)*f(rm(i),ym(j),zm(k));
                    end
                end
            end
            gz(jy,ix)=ss*rho*G*1e5;      
        end
    end

end