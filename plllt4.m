% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab code for plotting all tesseroid and Gaussfft based models and
%%comparing the results

clear all
close all

%all lon and lat intervals
dg_int=[0.8 1:15];
for idd=1:length(dg_int)
    %importing data

    %filename for tesseroid model
    flnm_data1=sprintf('gravity_fixed_density_tesseroid_%2.2f_degree.txt',dg_int(idd));
    data_tess=importdata(fullfile('.', 'output','tesseroid_data',flnm_data1));
    
    %filename for gaussfft model
    flnm_data2=sprintf('gravity_fixed_density_gauss_fft_%2.2f_degree.txt',dg_int(idd));
    data_gauss=importdata(fullfile('.', 'output','tesseroid_data',flnm_data2));
    
    %finding error 
    vv=abs(data_tess-data_gauss);
    max_error(idd,1)=max(vv(:)); min_error(idd,1)=min(vv(:)); 
    rel_rmse(idd,1)=(norm(vv)/norm(data_tess))*100; rel_ave_err(idd,1)=(mean(vv)/mean(abs(data_tess)))*100;
end
tbl=[dg_int' max_error min_error rel_rmse rel_ave_err];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%importing topography data for fixed density
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=linspace(72,72.8,129);
yy=linspace(8,8.8,129);

ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
ccmap2=makecolormap({'navy','lightseagreen','gold','magenta'}, 128);
ccmap3=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%plotting the topography
[XX,YY]=meshgrid(xx,yy);
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

xlabel('Lon')
ylabel('Lat')
zlabel('Depth(km)')
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

xlabel('Lon')
ylabel('Lat')
zlabel('Depth(km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])


%importing anomaly data for fixed density
gz=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_tesseroid_0.80_degree.txt')); 
xx1=linspace(72,72.8,128); yy1=linspace(8,8.8,128);
[XX1,YY1]=meshgrid(xx1,yy1);

%plotting the gravity anomalies
figure(1)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)

colormap(ax2,ccmap3)
shading interp
xlabel('Lon')
ylabel('Lat')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%importing topography data for fixed density
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=linspace(72,72+15,129);
yy=linspace(8,8+15,129);

ccmap1=makecolormap({'green','yellow','khaki','olive','brown','silver'}, 128);
ccmap2=makecolormap({'navy','lightseagreen','gold','magenta'}, 128);
ccmap3=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%plotting the topography
[XX,YY]=meshgrid(xx,yy);
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

xlabel('Lon')
ylabel('Lat')
zlabel('Depth(km)')
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

xlabel('Lon')
ylabel('Lat')
zlabel('Depth(km)')
%title('depth data plot')
set(gca, 'ZDir','reverse')
c = colorbar;
c.Label.String = 'Depth (km)';
xlim([min(XX(:)) max(XX(:))])
ylim([min(YY(:)) max(YY(:))])


%importing anomaly data for fixed density
gz=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_tesseroid_15.00_degree.txt')); 
xx1=linspace(72,72+15,128); yy1=linspace(8,8+15,128);
[XX1,YY1]=meshgrid(xx1,yy1);

%plotting the gravity anomalies
figure(2)
ax2=subplot(2,2,[1,3]);
surf(XX1,YY1,gz)

colormap(ax2,ccmap3)
shading interp
xlabel('Lon')
ylabel('Lat')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);


%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_tesseroid_0.80_degree.txt')); 
xx1=linspace(72,72.8,128); yy1=linspace(8,8.8,128);
[XX1,YY1]=meshgrid(xx1,yy1);

gz_gauss_fft=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_gauss_fft_0.80_degree.txt')); 

figure(3)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('Lon')
ylabel('Lat')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

dd1=abs(gz_gauss_fft-gz_true); 
%dd_min=0.5; dd_max=6;
figure(3)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd1)
%surf(XX1,YY1,abs(gz_analytical-gz_true))
shading interp
xlabel('Lon')
ylabel('Lat')
zlabel('Error (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);


%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_tesseroid_15.00_degree.txt')); 
xx1=linspace(72,72+15,128); yy1=linspace(8,8+15,128);
[XX1,YY1]=meshgrid(xx1,yy1);

gz_gauss_fft=importdata(fullfile('.', 'output','tesseroid_data','gravity_fixed_density_gauss_fft_15.00_degree.txt')); 

figure(4)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('Lon')
ylabel('Lat')
zlabel('Gravity anomaly (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

dd1=abs(gz_gauss_fft-gz_true); 
%dd_min=0.5; dd_max=6;
figure(4)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd1)
%surf(XX1,YY1,abs(gz_analytical-gz_true))
shading interp
xlabel('Lon')
ylabel('Lat')
zlabel('Error (mGal)')
set(gca,'TickDir','out');
view(2)
grid on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

