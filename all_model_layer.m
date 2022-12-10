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
%%surfaces and corresponding density contrats using vertically sliced
%%horizontal layers 
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
layers=150; %number of layers
tic
%gravity anomaly for given data
[XX1, YY1, gz]=grav_layer_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My,layers);
t1=toc;
fprintf('Computation time for exponential density layered model is %f\n',t1)

save(fullfile('.', 'output','gravity_exp_density_layer.txt'),'gz','-Ascii')
save(fullfile('.', 'output','x_meshgrid_exp_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_exp_density.txt'), 'YY1', '-Ascii')

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
layers=150; %number of layers
tic
%gravity anomaly for given data
[XX1, YY1, gz]=grav_layer_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My,layers);
t2=toc;
fprintf('Computation time for polynomial density layered model is %f\n',t2)

save(fullfile('.', 'output','gravity_polynomial_density_layer.txt'),'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_polynomial_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_polynomial_density.txt'), 'YY1', '-Ascii')

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
layers=150; %number of layers
tic
%gravity anomaly for given data
[XX1, YY1, gz]=grav_layer_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My,layers);
t3=toc;
fprintf('Computation time for complex density layered model is %f\n',t3)

save(fullfile('.', 'output','gravity_complex_density_layer.txt'),'gz', '-Ascii')
save(fullfile('.', 'output','x_meshgrid_complex_density.txt'), 'XX1', '-Ascii')
save(fullfile('.', 'output','y_meshgrid_complex_density.txt'), 'YY1', '-Ascii')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

