% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%Matlab comparision code for tesseroid model and gauss_fft based forward model for finding gravity 
%anomalies of topographic mass having fixed density distribution
clear all
close all
%all lon and lat intervals
dg_int=[0.8 1:30];
%loop for computing all tesseroid based gravity 
for idd=1:length(dg_int)
    %fixed density model
    %importing topography data
    data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
    data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));

    %Depth grids in degree
    lon=linspace(72,72+dg_int(idd),129); lat=linspace(8,8+dg_int(idd),129);
    m=length(lon); n=length(lat);
    [LON,LAT]=meshgrid(lon,lat);
    tic
    %observation point at z=0;
    Rad=0;
    %observation grids in meter
    
    %deeper depth
    [lon1,lat1,data2_g]=center_grid(lon,lat,data2);
    [LON1,LAT1]=meshgrid(lon1,lat1);

    %shallower depth
    [lon1,lat1,data1_g]=center_grid(lon,lat,data1);
    [LON1,LAT1]=meshgrid(lon1,lat1);

    rho=-400;  %density in kg/m^3 

    %number of Gaussian quadrature node
    N=2;
    %finding weights and node points
    [node,weight]=lgwt(N,-1,1);

   
%%
%Coputation for gaussfft models
    %lon lat and height from surface of Earth of topocentric coordinate
    lon=72; lat=8; h0=0; %lon lat in degree h0 in km)

    %span of the surface that has to plot in topocentric coordinate
    lon_span=linspace(lon,lon+dg_int(idd),129);
    lat_span=linspace(lat,lat+dg_int(idd),129);
    r_span=6371;

    %meshgrid of lon and lat
    [LON_SPAN,LAT_SPAN]=meshgrid(lon_span,lat_span);

    X_SPAN=r_span.*cos(deg2rad(LON_SPAN)).*cos(deg2rad(LAT_SPAN));
    Y_SPAN=r_span.*cos(deg2rad(LON_SPAN)).*sin(deg2rad(LAT_SPAN));
    Z_SPAN=r_span.*sin(deg2rad(LON_SPAN));

    %reference spheroid
    spheroid = referenceSphere('Earth');
    spheroid.LengthUnit = 'kilometer';

    sz=size(X_SPAN);
    %loop for finding x,y,z
    for i=1:sz(1)
        for j=1:sz(2)
            [yEast(i,j),xNorth(i,j),zUp(i,j)] = ecef2enu(X_SPAN(i,j),Y_SPAN(i,j),Z_SPAN(i,j),lon,lat,h0,spheroid);
        end
    end


    xx=(linspace(0,max(xNorth(:)),129)).*1000;
    yy=(linspace(0,max(yEast(:)),129)).*1000;

    %Gauss quadrature based estimation
    %observation point at z=0;
    z0=0;

    %density contrast
    rho=@(x,y,z) -400; 

    %number of gauss quadrature node
    Mx=4; My=4; 
    tic
    %gravity anomaly for given data
    [XX1, YY1, gz_gaussfft, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
    %[XX1, YY1, gz_stdfft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);
    t=toc;
    fprintf('For %2.2f degree interval, required computational time for gauss_fft is %f\n',dg_int(idd),t)
    flnm_data2=sprintf('gravity_fixed_density_gauss_fft_%2.2f_degree.txt',dg_int(idd));
    save(fullfile('.', 'output','tesseroid_data',flnm_data2),'gz_gaussfft', '-Ascii')
end

