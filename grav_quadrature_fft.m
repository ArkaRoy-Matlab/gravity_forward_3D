% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for finding gravity anomalies for topographic surfaces
%%having any 3D density contrast quadrature based standard FFT method

function [XX1, YY1, gz, delta1, delta2, N]=grav_quadrature_fft(data1,data2,xx,yy,rho,z0,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %input: 
    %   data1  = topography data for shallower surface in m
    %   data2  = topography data for deeper    surface in m
    %   xx     = 1D x grid locations in m
    %   yy     = 1D y grid locations in m
    %   rho    = density function in kg/m^3
    %   z0     = observation point in m
    %   L      = Grid expansion ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %output:
    %   XX1    = meshgrid of x observation points
    %   YY1    = meshgrid of y observation points
    %   gz     = gravity anomaly for given topgraphy and density contrasts
    %   delta1 = shallower surface mean depth for Taylor series expansion
    %   delta2 = deeper    surface mean depth for Taylor series expansion
    %   N      = number of Gaussian quadrature node for integral
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %prerequired functions
    %   1. lgwt.m
    %   2. find_delta.m
    %   3. find_nodes.m
    %   4. sfft2.m
    %   5. sifft2.m
    %   6. sfft_X.m
    %   7. sfft_Y.m
    %   8. sifft_X.m
    %   9. sittf_Y.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Depth grids in meter
    [XX,YY]=meshgrid(xx,yy);

    %observation grids in meter
    %deeper depth
    [xx1,yy1,data2_g]=center_grid(xx,yy,data2);
    [XX1,YY1]=meshgrid(xx1,yy1);

    %shallower depth
    [xx1,yy1,data1_g]=center_grid(xx,yy,data1);
    [XX1,YY1]=meshgrid(xx1,yy1);

    %data spacing along x and y direction
    dx=xx(2)-xx(1);
    dy=yy(2)-yy(1);

    %mean depth and data-mean_depth
    %deeper layer
    delta2=find_delta(data2_g,dx,dy);
    %shallower layer
    delta1=find_delta(data1_g,dx,dy);

    %%
    %gauss fft mass line model for calculating gravity anomalies for fixed density
    %Gauss quadrature node and weight
    coef1=(data1_g)./2;
    coef2=(data2_g)./2;
    %Gaussian quadrature node and weights for given density and topography
    N1=find_nodes(XX1,YY1,coef1,rho);
    N2=find_nodes(XX1,YY1,coef2,rho);
    N=max([N1 N2]);
    [y,w]=lgwt(N,-1,1);
    
    [m,n]=size(XX1); 
    %New length and width of data
    m_new=L*m;
    n_new=L*n;
    %frequency fourier transform
    kx=(2*pi/(m_new*dx)).*[(0:1:floor(m_new/2)-1) (floor(-m_new/2):1:-1)];
    ky=(2*pi/(n_new*dx)).*[(0:1:floor(n_new/2)-1) (floor(-n_new/2):1:-1)];


    %2d form of k
    [KX,KY]=meshgrid(kx,ky);
    %absolute value of k
    K=sqrt(KX.^2+KY.^2);
    %Change K with very small term to reduce divergence
    K(1,1)=10^-15;
    G=6.67408*10^-11; %in m^3kg^-1s^-2
    
    %gravity anomaly using Parker formula
    %%% forward gravity calculation
        hs1=2*pi*G.*exp(abs(K).*z0).*exp(-abs(K).*delta1);
        hs2=2*pi*G.*exp(abs(K).*z0).*exp(-abs(K).*delta2);
        tongF1=0; tongF2=0;
        for in=1:30
            nn=in-1;
            ss1=0; ss2=0;
            for jn=1:N
                ss1=ss1+w(jn).*rho(XX1,YY1,coef1*y(jn)+coef1).*(coef1*y(jn)+coef1-delta1).^nn;
                ss2=ss2+w(jn).*rho(XX1,YY1,coef2*y(jn)+coef2).*(coef2*y(jn)+coef2-delta2).^nn;
            end
            ss1=ss1.*coef1; ss2=ss2.*coef2;   
                   
            ss1_pad=padarray(ss1,[ceil(((L-1)/2)*n) ceil(((L-1)/2)*m)],0,'both');
            ss2_pad=padarray(ss2,[ceil(((L-1)/2)*n) ceil(((L-1)/2)*m)],0,'both');
            
            tongF1=tongF1+(((-abs(K)).^(nn))./(factorial(nn))).*fft2(ss1_pad);
            tongF2=tongF2+(((-abs(K)).^(nn))./(factorial(nn))).*fft2(ss2_pad);
        end
        
        Fg=hs2.*(tongF2).*10^5-hs1.*(tongF1).*10^5;
        g0=(ifft2(Fg));

    %%extracting data without padding
    gz=g0(ceil(((L-1)/2)*n)+1:ceil(((L-1)/2)*n)+n,ceil(((L-1)/2)*m)+1:ceil(((L-1)/2)*m)+m);
end