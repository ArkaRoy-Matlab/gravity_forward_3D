%Matlab code for gauss fft mass line model for finding gravity anomalies of topographic
%mass having fixed density distribution
clear all
close all

%fixed density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_fixed_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_fixed_density.txt'));

[XX,YY]=meshgrid(xx,yy);
tic
%observation point at z=0;
z0=0;
%observation grids in meter
%deeper depth
[xx1,yy1,data2_g]=center_grid(xx,yy,data2);
[XX1,YY1]=meshgrid(xx1,yy1);

%shallower depth
[xx1,yy1,data1_g]=center_grid(xx,yy,data1);
[XX1,YY1]=meshgrid(xx1,yy1);

%mean depth and data-mean_depth
%deeper layer
delta2=(max(data2_g(:))+min(data2_g(:)))/2;
deltah2=data2_g-delta2; 

%shallower layer
delta1=(max(data1_g(:))+min(data1_g(:)))/2;
deltah1=data1_g-delta1; 

%gauss fft mass line model for calculating gravity anomalies for fixed density
%data spacing along x and y direction
dx=xx(2)-xx(1);
dy=yy(2)-yy(1);
%Gauss quadrature node and weight
Mx=4; My=4;  %number of gauss quadrature node
[Abscissae_x,Wx]=lgwt(2*Mx,0,1);
[Abscissae_y,Wy]=lgwt(2*My,0,1);
n1=length(Wx);
n2=length(Wy)/2;
gz=0;
%gravitational constant
G=6.67408*10^-11; %in m^3kg^-1s^-2
%density distribution
rho=-400;

%dimension of new observation grid
m_new=length(xx1); n_new=length(yy1);
%wavevectors
dkx=2*pi/(m_new*dx);dky=2*pi/(n_new*dy);
nx1=ceil(-m_new/2);nx2=ceil(m_new/2-1);
ny1=ceil(-n_new/2);ny2=ceil(n_new/2-1);
NX=(nx1:1:nx2)';NY=(ny1:1:ny2)';% Column vector
[KX0,KY0]=meshgrid(dkx*NX,dky*NY);


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
        hs1=2*pi*G.*rho.*exp(abs(K).*z0).*exp(-abs(K).*delta1);
        hs2=2*pi*G.*rho.*exp(abs(K).*z0).*exp(-abs(K).*delta2);
        tongF1=0; tongF2=0;
        %Taylor series sum
        for im=1:30
              im=im-1;
              ss1=(deltah1).^(im); ss2=(deltah2).^(im);
              vv1=(sfft2(ss1./(factorial(im)),xi,eta,0.5,0.5));
              vv2=(sfft2(ss2./(factorial(im)),xi,eta,0.5,0.5));
               
             tongF1=tongF1+(((-abs(K)).^(im-1))).*vv1;
             tongF2=tongF2+(((-abs(K)).^(im-1))).*vv2;
        end
        Fg=hs2.*(tongF2).*10^5-hs1.*(tongF1).*10^5;
        temp_gz=sifft2(Fg,xi,eta,0.5,0.5);
        gz=gz+2*Wx(i1)*Wy(i2)*real(temp_gz);
    end   
end
t=toc;
%save the gravity anomay 
save(fullfile('.', 'output','gravity_fixed_density_analytic.txt'), 'gz', '-Ascii')