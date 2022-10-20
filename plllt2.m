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
%%corresponding gauss FFT and standard FFT quadrature based forward gravity anomalies
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ccmap1=makecolormap({'darkred','orange','khaki','thistle'}, 128);
ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);


%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','gravity_fixed_density_prism.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_fixed_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_fixed_density.txt')); 
XX1=XX1*10^-3; YY1=YY1*10^-3; 
gz_analytical=importdata(fullfile('.', 'output','gravity_fixed_density_analytic.txt')); 
gz_gauss_fft=importdata(fullfile('.', 'output','gravity_fixed_density_quadrature1.txt')); 
gz_fft=importdata(fullfile('.', 'output','gravity_fixed_density_quadrature2.txt')); 

figure(1)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_analytical)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

dd1=abs(gz_analytical-gz_true); 
dd2=abs(gz_gauss_fft-gz_true);
dd3=abs(gz_fft-gz_true);
dd_min=min([dd1(:);dd2(:);dd3(:)]); dd_max=max([dd1(:);dd2(:);dd3(:)]);
%dd_min=0.5; dd_max=6;
figure(1)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd1)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_analytical-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])


figure(2)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(2)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd2)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_gauss_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(3)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(3)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd3)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%All error calculations
%analytic
vv1=abs(gz_analytical-gz_true);
max_error_ana=max(vv1(:)); min_error_ana=min(vv1(:)); 
rel_rmse_ana=(norm(vv1)/norm(gz_true))*100; rel_ave_err_ana=(mean(vv1)/mean(abs(gz_true)))*100;

%Gauss-FFT quadrature
vv2=abs(gz_gauss_fft-gz_true);
max_error_gauss_fft=max(vv2(:)); min_error_gauss_fft=min(vv2(:)); 
rel_rmse_gauss_fft=(norm(vv2)/norm(gz_true))*100; rel_ave_err_gauss_fft=(mean(vv2)/mean(abs(gz_true)))*100;

%FFT quadrature
vv3=abs(gz_fft-gz_true);
max_error_FFT=max(vv3(:)); min_error_FFT=min(vv3(:)); 
rel_rmse_FFT=(norm(vv3)/norm(gz_true))*100; rel_ave_err_FFT=(mean(vv3)/mean(abs(gz_true)))*100;

fprintf('\nAll error calculations for model 1.\n')
fprintf('For analytic approach, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_ana,min_error_ana,rel_rmse_ana,rel_ave_err_ana)
fprintf('For Gauss-FFT quadrature, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_gauss_fft,min_error_gauss_fft,rel_rmse_gauss_fft,rel_ave_err_gauss_fft)
fprintf('For FFT quadrature, max error=%f, min error =%f , rel rmse error =%f and rel ave error =%f\n',max_error_FFT,min_error_FFT,rel_rmse_FFT,rel_ave_err_FFT)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%importing anomaly data for exponential density
gz_true=importdata(fullfile('.', 'output','gravity_exp_density_layer.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_exp_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_exp_density.txt')); 
XX1=XX1*10^-3; YY1=YY1*10^-3; 
gz_analytical=importdata(fullfile('.', 'output','gravity_exp_density_analytic.txt')); 
gz_gauss_fft=importdata(fullfile('.', 'output','gravity_exp_density_quadrature1.txt')); 
gz_fft=importdata(fullfile('.', 'output','gravity_exp_density_quadrature2.txt')); 

figure(4)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_analytical)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

dd1=abs(gz_analytical-gz_true); 
dd2=abs(gz_gauss_fft-gz_true);
dd3=abs(gz_fft-gz_true);
dd_min=min([dd1(:);dd2(:);dd3(:)]); dd_max=max([dd1(:);dd2(:);dd3(:)]);

figure(4)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd1)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_analytical-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])


figure(5)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(5)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd2)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_gauss_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(6)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(6)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd3)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%All error calculations
%analytic
vv1=abs(gz_analytical-gz_true);
max_error_ana=max(vv1(:)); min_error_ana=min(vv1(:)); 
rel_rmse_ana=(norm(vv1)/norm(gz_true))*100; rel_ave_err_ana=(mean(vv1)/mean(abs(gz_true)))*100;

%Gauss-FFT quadrature
vv2=abs(gz_gauss_fft-gz_true);
max_error_gauss_fft=max(vv2(:)); min_error_gauss_fft=min(vv2(:)); 
rel_rmse_gauss_fft=(norm(vv2)/norm(gz_true))*100; rel_ave_err_gauss_fft=(mean(vv2)/mean(abs(gz_true)))*100;

%FFT quadrature
vv3=abs(gz_fft-gz_true);
max_error_FFT=max(vv3(:)); min_error_FFT=min(vv3(:)); 
rel_rmse_FFT=(norm(vv3)/norm(gz_true))*100; rel_ave_err_FFT=(mean(vv3)/mean(abs(gz_true)))*100;

fprintf('\nAll error calculations for model 2.\n')
fprintf('For analytic approach, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_ana,min_error_ana,rel_rmse_ana,rel_ave_err_ana)
fprintf('For Gauss-FFT quadrature, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_gauss_fft,min_error_gauss_fft,rel_rmse_gauss_fft,rel_ave_err_gauss_fft)
fprintf('For FFT quadrature, max error=%f, min error =%f , rel rmse error =%f and rel ave error =%f\n',max_error_FFT,min_error_FFT,rel_rmse_FFT,rel_ave_err_FFT)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%importing anomaly data for polynomial density
gz_true=importdata(fullfile('.', 'output','gravity_polynomial_density_layer.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_polynomial_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_polynomial_density.txt')); 
XX1=XX1*10^-3; YY1=YY1*10^-3; 
gz_analytical=importdata(fullfile('.', 'output','gravity_polynomial_density_analytic.txt')); 
gz_gauss_fft=importdata(fullfile('.', 'output','gravity_polynomial_density_quadrature1.txt')); 
gz_fft=importdata(fullfile('.', 'output','gravity_polynomial_density_quadrature2.txt')); 
figure(7)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_analytical)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

dd1=abs(gz_analytical-gz_true); 
dd2=abs(gz_gauss_fft-gz_true);
dd3=abs(gz_fft-gz_true);
dd_min=min([dd1(:);dd2(:);dd3(:)]); dd_max=max([dd1(:);dd2(:);dd3(:)]);
dd_max=0.2;
figure(7)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd1)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits

shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])


figure(8)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(8)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd2)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits

shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(9)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(9)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd3)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%All error calculations
%analytic
vv1=abs(gz_analytical-gz_true);
max_error_ana=max(vv1(:)); min_error_ana=min(vv1(:)); 
rel_rmse_ana=(norm(vv1)/norm(gz_true))*100; rel_ave_err_ana=(mean(vv1)/mean(abs(gz_true)))*100;

%Gauss-FFT quadrature
vv2=abs(gz_gauss_fft-gz_true);
max_error_gauss_fft=max(vv2(:)); min_error_gauss_fft=min(vv2(:)); 
rel_rmse_gauss_fft=(norm(vv2)/norm(gz_true))*100; rel_ave_err_gauss_fft=(mean(vv2)/mean(abs(gz_true)))*100;

%FFT quadrature
vv3=abs(gz_fft-gz_true);
max_error_FFT=max(vv3(:)); min_error_FFT=min(vv3(:)); 
rel_rmse_FFT=(norm(vv3)/norm(gz_true))*100; rel_ave_err_FFT=(mean(vv3)/mean(abs(gz_true)))*100;

fprintf('\nAll error calculations for model 3.\n')
fprintf('For analytic approach, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_ana,min_error_ana,rel_rmse_ana,rel_ave_err_ana)
fprintf('For Gauss-FFT quadrature, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_gauss_fft,min_error_gauss_fft,rel_rmse_gauss_fft,rel_ave_err_gauss_fft)
fprintf('For FFT quadrature, max error=%f, min error =%f , rel rmse error =%f and rel ave error =%f\n',max_error_FFT,min_error_FFT,rel_rmse_FFT,rel_ave_err_FFT)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%importing anomaly data for complex density
gz_true=importdata(fullfile('.', 'output','gravity_complex_density_layer.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_complex_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_complex_density.txt')); 
XX1=XX1*10^-3; YY1=YY1*10^-3; 
gz_gauss_fft=importdata(fullfile('.', 'output','gravity_complex_density_quadrature1.txt')); 
gz_fft=importdata(fullfile('.', 'output','gravity_complex_density_quadrature2.txt')); 

dd2=abs(gz_gauss_fft-gz_true);
dd3=abs(gz_fft-gz_true);
dd_min=min([dd2(:);dd3(:)]); dd_max=max([dd2(:);dd3(:)]);

figure(10)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(10)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd2)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_gauss_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = 'log_{10}(Error) (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(11)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

figure(11)
ax2=subplot(1,2,2);
surf(XX1,YY1,dd3)
caxis(ax2,[dd_min,dd_max]);    % set colorbar limits
%surf(XX1,YY1,abs(gz_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = 'log_{10}(Error) (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%All error calculations


%Gauss-FFT quadrature
vv2=abs(gz_gauss_fft-gz_true);
max_error_gauss_fft=max(vv2(:)); min_error_gauss_fft=min(vv2(:)); 
rel_rmse_gauss_fft=(norm(vv2)/norm(gz_true))*100; rel_ave_err_gauss_fft=(mean(vv2)/mean(abs(gz_true)))*100;

%FFT quadrature
vv3=abs(gz_fft-gz_true);
max_error_FFT=max(vv3(:)); min_error_FFT=min(vv3(:)); 
rel_rmse_FFT=(norm(vv3)/norm(gz_true))*100; rel_ave_err_FFT=(mean(vv3)/mean(abs(gz_true)))*100;

fprintf('\nAll error calculations for model 4.\n')
fprintf('For Gauss-FFT quadrature, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_gauss_fft,min_error_gauss_fft,rel_rmse_gauss_fft,rel_ave_err_gauss_fft)
fprintf('For FFT quadrature, max error=%f, min error =%f , rel rmse error =%f and rel ave error =%f\n',max_error_FFT,min_error_FFT,rel_rmse_FFT,rel_ave_err_FFT)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
ccmap1=makecolormap({'thistle','khaki','orange','darkred'}, 128);
ccmap2=makecolormap({'navy','firebrick','red','orange','violet','darkorchid','mediumspringgreen'}, 128);
%importing anomaly data for real sedimentary basin
gz_true=importdata(fullfile('.', 'output','gravity_real_santos.txt')); 
XX1=importdata(fullfile('.', 'output','x_meshgrid_real_density.txt')); 
YY1=importdata(fullfile('.', 'output','y_meshgrid_real_density.txt')); 
XX1=XX1*10^-3; YY1=YY1*10^-3; 
gz_gauss_fft=importdata(fullfile('.', 'output','gravity_real_density_quadrature.txt')); 
figure(12)
ax1=subplot(1,2,1);
surf(XX1,YY1,gz_gauss_fft)
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Gravity anomaly (mGal)')
view(2)
grid on;
box on;
colormap(ax1,ccmap2)
c = colorbar;
c.Label.String = 'Gravity anomaly (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])
figure(12)
ax2=subplot(1,2,2);
surf(XX1,YY1,abs(gz_gauss_fft-gz_true))
shading interp
xlabel('x (km)')
ylabel('y (km)')
zlabel('Error (mGal)')
view(2)
grid on;
box on;
colormap(ax2,ccmap1)
c = colorbar;
c.Label.String = '|Error| (mGal)';
xlim([min(XX1(:)) max(XX1(:))])
ylim([min(YY1(:)) max(YY1(:))])

%All error calculations
%Gauss-FFT quadrature
vv2=abs(gz_gauss_fft-gz_true);
max_error_gauss_fft=max(vv2(:)); min_error_gauss_fft=min(vv2(:)); 
rel_rmse_gauss_fft=(norm(vv2)/norm(gz_true))*100; rel_ave_err_gauss_fft=(mean(vv2)/mean(abs(gz_true)))*100;

fprintf('\nAll error calculations for real model.\n')
fprintf('For Gauss-FFT quadrature, max error=%f, min error =%2.2e , rel rmse error =%f and rel ave error =%f\n',max_error_gauss_fft,min_error_gauss_fft,rel_rmse_gauss_fft,rel_ave_err_gauss_fft)
