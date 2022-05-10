% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for finding gravity anomalies for topographic surfaces
%%having any 3D density contrast using vertially sliced horizontal layers

function [XX1, YY1, g_sum]=grav_layer_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My,nl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %input: 
    %   data1  = topography data for shallower surface in m
    %   data2  = topography data for deeper    surface in m
    %   xx     = 1D x grid locations in m
    %   yy     = 1D y grid locations in m
    %   rho    = density function in kg/m^3
    %   z0     = observation point in m
    %   Mx     = nodes for gauss fft in x direction (default is 2)
    %   My     = nodes for gauss fft in y direction (default is 2)
    %   nl     = number of layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %output:
    %   XX1    = meshgrid of x observation points
    %   YY1    = meshgrid of y observation points
    %   gz     = gravity anomaly for given topgraphy and density contrasts
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %prerequired functions
    %   1. lgwt.m
    %   2. sfft2.m
    %   3. sifft2.m
    %   4. sfft_X.m
    %   5. sfft_Y.m
    %   6. sifft_X.m
    %   7. sittf_Y.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Depth grids in meter
    [XX,YY]=meshgrid(xx,yy);

    %observation grids in meter
%deeper depth
[xx1,yy1,data2_g]=center_grid(xx,yy,data2);
[XX1,YY1]=meshgrid(xx1,yy1);

%shallower depth
[xx1,yy1,data1_g]=center_grid(xx,yy,data1);
[XX1,YY1]=meshgrid(xx1,yy1);

%gauss fft mass line model for calculating gravity anomalies for fixed density
%data spacing along x and y direction
dx=xx(2)-xx(1);
dy=yy(2)-yy(1);
%Gauss quadrature node and weight
[Abscissae_x,Wx]=lgwt(2*Mx,0,1);
[Abscissae_y,Wy]=lgwt(2*My,0,1);
n1=length(Wx);
n2=length(Wy)/2;
gz=0;
%gravitational constant
G=6.67408*10^-11; %in m^3kg^-1s^-2

%dimension of new observation grid
m_new=length(xx1); n_new=length(yy1);
%wavevectors
dkx=2*pi/(m_new*dx);dky=2*pi/(n_new*dy);
nx1=ceil(-m_new/2);nx2=ceil(m_new/2-1);
ny1=ceil(-n_new/2);ny2=ceil(n_new/2-1);
NX=(nx1:1:nx2)';NY=(ny1:1:ny2)';% Column vector
[KX0,KY0]=meshgrid(dkx*NX,dky*NY);

%%all depth layers 
layers=linspace(0,max(data2_g(:)),nl);
g_sum=0;
%loop for finding gravity anomalies as per depth layer bounds
for ll=1:length(layers)-1
    ll
    z1=layers(ll);  %shallower layer
    z2=layers(ll+1);%deeper layer
    zz1=data1_g;
    zz2=data2_g;
    
    zz1(zz1>z2)=z2; zz1(zz1<z1)=z1;
    zz2(zz2>z2)=z2; zz2(zz2<z1)=z1;
    
    delta1=mean(zz1(:)); %mean depth of shallower layer
    deltah1=zz1-delta1;
    
    delta2=mean(zz2(:)); %mean depth
    deltah2=zz2-delta2;
    %max(zz2(:))/delta2
    mn_depth=(z1+z2)/2;
    gz=0;
    %loop for gauss quadrature integral
    for i2=1:1:n2
        eta=Abscissae_y(i2);
        wy=Wy(i2);
        for i1=1:1:n1
            xi=Abscissae_x(i1);
            wx=Wx(i1);
            KX=KX0+xi*dkx;
            KY=KY0+eta*dky;
            K=sqrt(KX.^2+KY.^2);
            %gravity anomaly using Parker formula
            %%% forward gravity calculation
            hs1=2*pi*G.*exp(abs(K).*z0).*exp(-abs(K).*delta1);
            hs2=2*pi*G.*exp(abs(K).*z0).*exp(-abs(K).*delta2);
            tongF1=0; tongF2=0;
            %Taylor series sum
            for im=1:20
                 im=im-1;
                  ss1=squeeze(rho(XX1,YY1,mn_depth)).*(deltah1).^(im); ss2=squeeze(rho(XX1,YY1,mn_depth)).*(deltah2).^(im); 
                  vv1=(sfft2(ss1,xi,eta,0.5,0.5));
                  vv2=(sfft2(ss2,xi,eta,0.5,0.5));

                 tongF1=tongF1+(((-abs(K)).^(im-1))./(factorial(im))).*vv1;
                 tongF2=tongF2+(((-abs(K)).^(im-1))./(factorial(im))).*vv2;
            end
            Fg=hs2.*(tongF2).*10^5-hs1.*(tongF1).*10^5;
            temp_gz=sifft2(Fg,xi,eta,0.5,0.5);
            gz=gz+2*Wx(i1)*Wy(i2)*(temp_gz);
        end

    end
    gz=(real(gz));
    g_sum=g_sum+gz;
end
