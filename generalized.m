%%A generalized Matlab code for finding gravity anomalies for different topographic
%%surfaces and corresponding density contrats using Gauss-FFT and standard
%%FFT quadrature based model and layered model
clear all
close all

%%%Inputs%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exponential density model
%importing topography data

%shallower surface
data1=importdata(fullfile('.', 'input','synthetic_topo_exp_density_shallower_layer.txt'));
%deeper surface
data2=importdata(fullfile('.', 'input','synthetic_topo_exp_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_exp_density.txt')); %along x axis
yy=importdata(fullfile('.', 'input','synthetic_y_exp_density.txt')); %along y axis

%observation point at z=0;
z0=0;

%density contrast
rho=@(x,y,z) -500.*(2.32.*10^-5.*x+1.5.*10^-5.*y).*exp(-0.0187.*z.*10^-2);     %exponential

%number of gauss quadrature node
Mx=2; My=2; 

%number of layers for grav_layer_gaussfft
layers=150; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plotting the input data
%colormaps
ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(1)

surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
surf(XX,YY,data2./10^3)
[C,h] = contour3(XX,YY,data2./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ccmap1)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
grid on;
box on;
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Ouputs%%%
%%Gauss-fft quadrature based model%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gravity anomaly for given data using gauss-fft quadrature based model
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);

%saving the output data
save(fullfile('.', 'output','gravity_exp_density_quadrature1.txt'), 'gz', '-Ascii') %gravity anomaly
save(fullfile('.', 'output','x_meshgrid_exp_density.txt'), 'XX1', '-Ascii') %x grid
save(fullfile('.', 'output','y_meshgrid_exp_density.txt'), 'YY1', '-Ascii') %y grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting the output data
%colormaps
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);

XX1=XX1*10^-3; YY1=YY1*10^-3; 
figure(2)
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Standard fft quadrature based model%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gravity anomaly for given data using standard fft quadrature based model
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0); %gravity anomaly

%saving the output data
save(fullfile('.', 'output','gravity_exp_density_quadrature2.txt'),'gz','-Ascii')
save(fullfile('.', 'output','x_meshgrid_exp_density.txt'), 'XX1', '-Ascii') %x grid
save(fullfile('.', 'output','y_meshgrid_exp_density.txt'), 'YY1', '-Ascii') %y grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting the output data
%colormaps
ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);

XX1=XX1*10^-3; YY1=YY1*10^-3; 
figure(3)
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%layer based model%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XX1, YY1, gz]=grav_layer_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My,layers); 
%gravity anomaly for given data using layer based model
%saving the output data
save(fullfile('.', 'output','gravity_exp_density_layer.txt'),'gz','-Ascii')%gravity anomaly
save(fullfile('.', 'output','x_meshgrid_exp_density.txt'), 'XX1', '-Ascii')%x grid
save(fullfile('.', 'output','y_meshgrid_exp_density.txt'), 'YY1', '-Ascii')%y grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%colormaps
ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);

XX1=XX1*10^-3; YY1=YY1*10^-3; 
figure(4)
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
