%%Matlab code for coordinate transform 
clear all
close all

%lon lat and height from surface of Earth of topocentric coordinate
lon=72; lat=8; h0=0; %lon lat in degree h0 in km)

%span of the surface that has to plot in topocentric coordinate
lon_span=linspace(72,72.8,129);
lat_span=linspace(8,8.8,129);
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

surf(xNorth,yEast,zUp)
%surf(LON_SPAN,LAT_SPAN,zUp)

xx=linspace(0,max(xNorth(:)),129);
yy=linspace(0,max(yEast(:)),129);

[X,Y]=meshgrid(xx,yy);
Z_mesh=griddata(xNorth,yEast,(zUp),X,Y);
Z_mesh(isnan(Z_mesh))=0;
figure(2)
surf(X,Y,Z_mesh)

%fixed density model
%importing topography data
data1=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_shallower_layer.txt'));
data2=importdata(fullfile('.', 'input','synthetic_topo_fixed_density_deeper_layer.txt'));
data1=data1./15; data2=data2./15;
[x_g,y_g,data_Z]=center_grid(xx,yy,abs(Z_mesh));
%observation point at z=0;
%z0=data_Z*1000;
z0=0; 
%density contrast
rho=@(x,y,z) -400; 
xx=xx.*1000; yy=yy.*1000;
%xx=(lon-lon(1)).*111.321*1000; yy=(lat-lat(1)).*111.321*1000;
%number of gauss quadrature node
Mx=4; My=4; 
tic
%gravity anomaly for given data
[XX1, YY1, gz_gaussfft, delta1, delta2, N]=grav_quadrature_gaussfft(data1,data2,xx,yy,rho,z0,Mx,My);
%[XX1, YY1, gz_stdfft, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,8);

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
    ii
    for jj=1:n-1
        
        xp1=xx(ii)*10^-3; xp2=xx(ii+1)*10^-3;  %x grid for each prism
        yp1=yy(jj)*10^-3; yp2=yy(jj+1)*10^-3;  %y gris for each prism
        %mean depth for each prism 
        zp1=data1_g(ii,jj)*10^-3; %shallower depth
        zp2=data2_g(ii,jj)*10^-3; %deeper depth
        gz=gz+gprism(xp1,yp1,zp1,xp2,yp2,zp2,xx1*10^-3,yy1*10^-3,z0,rho);
    end
end
gz=gz';

save(fullfile('.', 'output','gravity_fixed_density_gaussfft_comparision.txt'),'gz', '-Ascii')
%save(fullfile('.', 'output','gravity_fixed_density_stdfft_comparision.txt'),'gz_stdfft', '-Ascii')

save(fullfile('.', 'output','x_meshgrid_fixed_density_comparision.txt'),'XX1','-Ascii')


figure(3)
surf(gz)
shading interp

gz_tess=importdata(fullfile('.', 'output','gravity_fixed_density_tesseroid_comparision.txt'));
figure(4)
surf(gz_tess)
shading interp
figure(5)
surf(gz_tess-gz)

figure(6)
surf(data1)
hold on 
surf(data2)
shading interp