% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%Matlab code for prism model for finding gravity anomalies of topographic
%mass having fixed density distribution using prismatic approach
clear all
close all

%fixed density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_fixed_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_fixed_density.txt'));
m=length(xx); n=length(yy);
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

rho=-400;  %density in kg/m^3 
gz=0;
%loop for finding gravity anomalies for each prism
for ii=1:m-1
    for jj=1:n-1
        
        xp1=xx(ii)*10^-3; xp2=xx(ii+1)*10^-3;  %x grid for each prism
        yp1=yy(jj)*10^-3; yp2=yy(jj+1)*10^-3;  %y gris for each prism
        %mean depth for each prism 
        zp1=data1_g(jj,ii)*10^-3; %shallower depth
        zp2=data2_g(jj,ii)*10^-3; %deeper depth
        gz=gz+gprism(xp1,yp1,zp1,xp2,yp2,zp2,xx1*10^-3,yy1*10^-3,z0,rho);
    end
end
t=toc;
fprintf('Computation time for fixed density prism model is %f\n',t)
save(fullfile('.', 'output','gravity_fixed_density_prism.txt'),'gz', '-Ascii')