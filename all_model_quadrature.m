% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab code for finding gravity anomalies for different topographic
%%surfaces and corresponding density contrats using Gauss-FFT and standard
%%FFT quadrature based model
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_fixed_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_fixed_density.txt'));

%observation point at z=0;
z0=0;

%density contrast
rho=@(x,y,z) -400; 

%number of gauss quadrature node
Mx=2; My=2; 
tic
%gravity anomaly for given data
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
t11=toc;
fprintf('For model 1 gauss fft computational time =%f\n',t11)
fprintf('For fixed density model\n')
fprintf('\t delta1=%f m and delta2=%f m\n',delta1,delta2)
fprintf('\t Gaussian quadrature integral nodes =%d\n',N)
save(fullfile('.', 'output','gravity_fixed_density_quadrature1.txt'),'gz','-Ascii')
save(fullfile('.', 'output','x_meshgrid_fixed_density.txt'),'XX1','-Ascii')
save(fullfile('.', 'output','y_meshgrid_fixed_density.txt'),'YY1','-Ascii')
tic
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);
t12=toc;
save(fullfile('.', 'output','gravity_fixed_density_quadrature2.txt'),'gz','-Ascii')
fprintf('For model 1 standard fft computational time =%f\n',t12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%exponential density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_exp_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_exp_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_exp_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_exp_density.txt'));

%observation point at z=0;
z0=0;

%density contrast
rho=@(x,y,z) -500.*(2.32.*10^-5.*x+1.5.*10^-5.*y).*exp(-0.0187.*z.*10^-2);     %exponential

%number of gauss quadrature node
Mx=2; My=2; 
tic
%gravity anomaly for given data
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
t21=toc;
fprintf('For model 2 gauss fft computational time =%f\n',t21)
fprintf('For exponential density model\n')
fprintf('\t delta1=%f m and delta2=%f m\n',delta1,delta2)
fprintf('\t Gaussian quadrature integral nodes =%d\n',N)

save(fullfile('.', 'output','gravity_exp_density_quadrature1.txt'), 'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_exp_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_exp_density.txt'), 'YY1', '-Ascii')
tic
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);
t22=toc;
save(fullfile('.', 'output','gravity_exp_density_quadrature2.txt'),'gz','-Ascii')
fprintf('For model 2 standard fft computational time =%f\n',t22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%polynomial density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_polynomial_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_polynomial_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_polynomial_density.txt'));

%observation point at z=0;
z0=0;

%density contrast
%rho=@(x,y,z) -300-0.3435.*10^-5.*z-0.6764.*10^-7.*z.^2-0.04247.*10^-11.*z.^3;  %polynomial
rho=@(x,y,z) -300-0.1435.*10^-3.*z-0.1764.*10^-5.*z.^2-0.4247.*10^-10.*z.^3;  %polynomial
%number of gauss quadrature node
Mx=2; My=2; 
tic
%gravity anomaly for given data
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
t31=toc;
fprintf('For model 3 gauss fft computational time =%f\n',t31)
fprintf('For polynomial density model\n')
fprintf('\t delta1=%f m and delta2=%f m\n',delta1,delta2)
fprintf('\t Gaussian quadrature integral nodes =%d\n',N)

save(fullfile('.', 'output','gravity_polynomial_density_quadrature1.txt'), 'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_polynomial_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_polynomial_density.txt'), 'YY1', '-Ascii')
tic
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);
t32=toc;
save(fullfile('.', 'output','gravity_polynomial_density_quadrature2.txt'), 'gz', '-Ascii')
fprintf('For model 3 standard fft computational time =%f\n',t32)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%complex density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_complex_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_complex_density_deeper_layer.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','synthetic_x_complex_density.txt'));
yy=importdata(fullfile('.', 'input','synthetic_y_complex_density.txt'));

%observation point at z=0;
z0=0;

%density contrast
rho=@(x,y,z) ((7^3)./(-0.25-0.1711.*z.^2.*10^-7))+(163+6.36*10^-7*x).*cos(3.2+9*10^-7.*y);  %hyperbolic

%number of gauss quadrature node
Mx=2; My=2; 
tic
%gravity anomaly for given data
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
t41=toc;
fprintf('For model 4 gauss fft computational time =%f\n',t41)
fprintf('Forcomplex density model\n')
fprintf('\t delta1=%f m and delta2=%f m\n',delta1,delta2)
fprintf('\t Gaussian quadrature integral nodes =%d\n',N)

save(fullfile('.', 'output','gravity_complex_density_quadrature1.txt'), 'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_complex_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_complex_density.txt'), 'YY1', '-Ascii')
tic
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);
t42=toc;
save(fullfile('.', 'output','gravity_complex_density_quadrature2.txt'),'gz','-Ascii')
fprintf('For model 4 standard fft computational time =%f\n',t42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Real density model
%importing topography data
data1=importdata(fullfile('.', 'input','real_topo_santos.txt'));
data1=0.*data1;
data2=importdata(fullfile('.', 'input','real_topo_santos.txt'));

%Depth grids in meter
xx=importdata(fullfile('.', 'input','real_x_santos.txt'));
yy=importdata(fullfile('.', 'input','real_y_santos.txt'));

%observation point at z=0;
z0=0;

%density contrast
rho=@(x,y,z) ((-0.5.^3)./((-0.5-0.1711.*z.*10^-3).^2)).*1000;  %parabolic

%number of gauss quadrature node
Mx=2; My=2; 
tic
%gravity anomaly for given data
[XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
t51=toc;
fprintf('For real model gauss fft computational time =%f\n',t51)
fprintf('For real density model\n')
fprintf('\t delta1=%f m and delta2=%f m\n',delta1,delta2)
fprintf('\t Gaussian quadrature integral nodes =%d\n',N)

save(fullfile('.', 'output','gravity_real_density_quadrature.txt'), 'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_real_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_real_density.txt'), 'YY1', '-Ascii')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


