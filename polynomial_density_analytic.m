%Matlab code for gauss fft mass line model for finding gravity anomalies of topographic
%mass having polynomial density distribution
clear all
close all

%fixed density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_polynomial_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_polynomial_density.txt'));

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
%rho=@(x,y,z) -300-0.3435.*10^-5.*z-0.6764.*10^-7.*z.^2-0.04247.*10^-11.*z.^3;  %polynomial
rho=@(x,y,z) -300-0.1435.*10^-3.*z-0.1764.*10^-5.*z.^2-0.4247.*10^-10.*z.^3;  %polynomial
%dimension of new observation grid
m_new=length(xx1); n_new=length(yy1);
%wavevectors
dkx=2*pi/(m_new*dx);dky=2*pi/(n_new*dy);
nx1=ceil(-m_new/2);nx2=ceil(m_new/2-1);
ny1=ceil(-n_new/2);ny2=ceil(n_new/2-1);
NX=(nx1:1:nx2)';NY=(ny1:1:ny2)';% Column vector
[KX0,KY0]=meshgrid(dkx*NX,dky*NY);
%a0=-300; a1=-0.3435.*10^-5; 
%a2=-0.6764.*10^-7; a3=-0.04247.*10^-11;

a0=-300; a1=-0.1435.*10^-3; 
a2=-0.1764.*10^-5; a3=-0.4247.*10^-10;

A1=(a0+a1.*delta1+a2.*delta1.^2+a3.*delta1.^3);
A2=(a1+2*a2*delta1+3*a3*delta1^2);
A3=(a2+3*a3*delta1);
A4=a3;

B1=(a0+a1.*delta2+a2.*delta2.^2+a3.*delta2.^3);
B2=(a1+2*a2*delta2+3*a3*delta2^2);
B3=(a2+3*a3*delta2);
B4=a3;

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
        for im=1:36
              im=im-1;
              ss11=A1.*((deltah1).^(im+1)-(-delta1).^(im+1))./(im+1); ss21=B1.*((deltah2).^(im+1)-(-delta2).^(im+1))./(im+1);
              ss12=A2.*((deltah1).^(im+2)-(-delta1).^(im+2))./(im+2); ss22=B2.*((deltah2).^(im+2)-(-delta2).^(im+2))./(im+2);
              ss13=A3.*((deltah1).^(im+3)-(-delta1).^(im+3))./(im+3); ss23=B3.*((deltah2).^(im+3)-(-delta2).^(im+3))./(im+3);
              ss14=A4.*((deltah1).^(im+4)-(-delta1).^(im+4))./(im+4); ss24=B4.*((deltah2).^(im+4)-(-delta2).^(im+4))./(im+4);
              %ss1=-(-delta1)^im; ss2=-(-delta2)^im;
              vv1=sfft2((ss11+ss12+ss13+ss14),xi,eta,0.5,0.5);
              vv2=sfft2((ss21+ss22+ss23+ss24),xi,eta,0.5,0.5);
               
             tongF1=tongF1+(((-abs(K)).^(im))).*vv1./(factorial(im));
             tongF2=tongF2+(((-abs(K)).^(im))).*vv2./(factorial(im));
             
             
        end
        Fg=hs2.*(tongF2).*10^5-hs1.*(tongF1).*10^5;
        temp_gz=sifft2(Fg,xi,eta,0.5,0.5);
        gz=gz+2*Wx(i1)*Wy(i2)*real(temp_gz);
    end   
end
t=toc;
fprintf('Computation time for polynomial density analytic model is %f\n',t)
%save the gravity anomay 
save(fullfile('.', 'output','gravity_polynomial_density_analytic.txt'), 'gz', '-Ascii')