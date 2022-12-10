% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab code for plotting all synthetic topographic surfaces and
%%corresponding analytic forward gravity anomalies
clear all
close all
%importing topography data for fixed density
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_fixed_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_fixed_density.txt'));

ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
ccmap2=makecolormap({'navy','lightseagreen','gold','magenta'}, 128);
ccmap3=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%plotting the topography
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(1)
ax1=subplot(2,2,2);
%plotting the data
surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap1)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

ax1=subplot(2,2,4);
%plotting the data
surf(XX,YY,data2./10^3)
hold on
[C,h] = contour3(XX,YY,data2./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap2)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])


%importing anomaly data for fixed density
gz=importdata(fullfile('.', 'output','gravity_fixed_density_prism.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_fixed_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_fixed_density.txt')); 

XX1=XX1./1000; YY1=YY1./1000;

%plotting the gravity anomalies
figure(1)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)

colormap(ax2,ccmap3)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
%
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%title('gravity anomaly gauss fft model')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%importing topography data for exponential density
data1=importdata(fullfile('.', 'input','synthetic_topo_exp_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_exp_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_exp_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_exp_density.txt'));

%plotting the topography
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(2)
ax1=subplot(2,2,2);
%plotting the data
surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap1)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

ax1=subplot(2,2,4);
%plotting the data
surf(XX,YY,data2./10^3)
hold on
[C,h] = contour3(XX,YY,data2./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap2)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])


%importing anomaly data for fixed density
gz=importdata(fullfile('.', 'output','gravity_exp_density_analytic.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_exp_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_exp_density.txt')); 

XX1=XX1./1000; YY1=YY1./1000;

%plotting the gravity anomalies
figure(2)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)

colormap(ax2,ccmap3)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
%
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%importing topography data for polynomial density
data1=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_polynomial_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_polynomial_density.txt'));

%plotting the topography
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(3)
ax1=subplot(2,2,2);
%plotting the data
surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap1)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

ax1=subplot(2,2,4);
%plotting the data
surf(XX,YY,data2./10^3)
hold on
[C,h] = contour3(XX,YY,data2./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap2)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])


%importing anomaly data
gz=importdata(fullfile('.', 'output','gravity_polynomial_density_layer.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_polynomial_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_polynomial_density.txt')); 

XX1=XX1./1000; YY1=YY1./1000;

%plotting the gravity anomalies
figure(3)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax2,ccmap3)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])
%title('gravity anomaly gauss fft model')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%importing topography data for complex density
data1=importdata(fullfile('.', 'input','synthetic_topo_complex_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_complex_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_complex_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_complex_density.txt'));

%plotting the topography
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(4)
ax1=subplot(2,2,2);
%plotting the data
surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap1)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

ax1=subplot(2,2,4);
%plotting the data
surf(XX,YY,data2./10^3)
hold on
[C,h] = contour3(XX,YY,data2./10^3,20,'k');
shading interp
set(gca, 'Zdir', 'reverse')

colormap(ax1,ccmap2)
view(3);

xlabel('x (km)')
ylabel('y (km)')
zlabel('Depth (km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])

%importing anomaly data
gz=importdata(fullfile('.', 'output','gravity_complex_density_layer.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_complex_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_complex_density.txt')); 

XX1=XX1./1000; YY1=YY1./1000;

%plotting the gravity anomalies
figure(4)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax2,ccmap3)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])
%title('gravity anomaly gauss fft model')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%importing topography data for real sedimentary basin
data1=importdata(fullfile('.', 'input','real_topo_santos.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','real_x_santos.txt'));
yy=importdata(fullfile('.', 'input','real_y_santos.txt'));

%plotting the topography
[XX,YY]=meshgrid(xx,yy);
XX=XX./1000; YY=YY./1000;
figure(5)
ax1=subplot(2,1,2);
%plotting the data
surf(XX,YY,data1./10^3)
hold on;
[C,h] = contour3(XX,YY,data1./10^3,20,'k');
surf(XX,YY,data1./10^3)
shading interp
set(gca, 'Zdir', 'reverse')
colormap(ax1,ccmap1)
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

%importing anomaly data
gz=importdata(fullfile('.', 'output','gravity_real_santos.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_real_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_real_density.txt')); 

XX1=XX1./1000; YY1=YY1./1000;

%plotting the gravity anomalies
figure(5)
ax2=subplot(2,1,1);
surf(XX1,YY1,gz)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax2,ccmap3)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])
%title('gravity anomaly gauss fft model')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

