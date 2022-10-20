% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab code for sensitivity study for finding gravity anomalies using standard
%%FFT quadrature based model for different grid expansion ratios
clear all
close all

%importing all data 
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

%True anomalies
%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','gravity_fixed_density_prism.txt')); 

%loop for different grid expansion ratios
for L=1:15
    %gravity anomaly for given data
    tic
    [XX1, YY1, gz_fft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,L);
    t(L,1)=toc;
    %FFT quadrature
    vv=abs(gz_fft-gz_true);
    max_error_FFT(L,1)=max(vv(:)); min_error_FFT(L,1)=min(vv(:)); 
    rel_rmse_FFT(L,1)=(norm(vv)/norm(gz_true))*100; rel_ave_err_FFT(L,1)=(mean(vv)/mean(abs(gz_true)))*100;
    fprintf('Model1 %d iteration completed\n',L)
end
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

%True anomalies
%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','gravity_exp_density_layer.txt')); 

%loop for different grid expansion ratios
for L=1:15
    %gravity anomaly for given data
    tic
    [XX1, YY1, gz_fft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,L);
    t(L,2)=toc;
    %FFT quadrature
    vv=abs(gz_fft-gz_true);
    max_error_FFT(L,2)=max(vv(:)); min_error_FFT(L,2)=min(vv(:)); 
    rel_rmse_FFT(L,2)=(norm(vv)/norm(gz_true))*100; rel_ave_err_FFT(L,2)=(mean(vv)/mean(abs(gz_true)))*100;
    fprintf('Model1 %d iteration completed\n',L)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
rho=@(x,y,z) -300-0.3435.*10^-5.*z-0.6764.*10^-7.*z.^2-0.04247.*10^-11.*z.^3;  %polynomial

%True anomalies
%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','gravity_polynomial_density_layer.txt')); 

%loop for different grid expansion ratios
for L=1:15
    %gravity anomaly for given data
    tic
    [XX1, YY1, gz_fft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,L);
    t(L,3)=toc;
    %FFT quadrature
    vv=abs(gz_fft-gz_true);
    max_error_FFT(L,3)=max(vv(:)); min_error_FFT(L,3)=min(vv(:)); 
    rel_rmse_FFT(L,3)=(norm(vv)/norm(gz_true))*100; rel_ave_err_FFT(L,3)=(mean(vv)/mean(abs(gz_true)))*100;
    fprintf('Model1 %d iteration completed\n',L)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
rho=@(x,y,z) ((18.5^3)./(-0.25-0.1711.*z.^2.*10^-7))+(163+6.36*10^-7*x).*cos(3.2+9*10^-7.*y);  %hyperbolic

%True anomalies
%importing anomaly data for fixed density
gz_true=importdata(fullfile('.', 'output','gravity_complex_density_layer.txt')); 

%loop for different grid expansion ratios
for L=1:15
    %gravity anomaly for given data
    tic
    [XX1, YY1, gz_fft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,L);
    t(L,4)=toc;
    %FFT quadrature
    vv=abs(gz_fft-gz_true);
    max_error_FFT(L,4)=max(vv(:)); min_error_FFT(L,4)=min(vv(:)); 
    rel_rmse_FFT(L,4)=(norm(vv)/norm(gz_true))*100; rel_ave_err_FFT(L,4)=(mean(vv)/mean(abs(gz_true)))*100;
    fprintf('Model1 %d iteration completed\n',L)
end

save(fullfile('.', 'output','sensitivity_time.txt'), 't', '-Ascii')
save(fullfile('.', 'output','sensitivity_rmse.txt'), 'rel_rmse_FFT', '-Ascii')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%