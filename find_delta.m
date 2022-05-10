% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for finding best delta for any topographic surface for
%%approximating exponential term in Taylor series used in Gauss-FFT forward
%%Modelling

function delta=find_delta(data,dx,dy)
%input:
    %   data=gridded data for topographic surface
    %   dx=grid spacing in x direction
    %   dy=grid spacing in y direction
    
    %output
    %   delta=delta value for approximating exponential term in Taylor series
    
    %data size
    sz=size(data);
    m_new=sz(1); n_new=sz(2);
    %wavevectors
    dkx=2*pi/(m_new*dx);dky=2*pi/(n_new*dy);
    nx1=ceil(-m_new/2);nx2=ceil(m_new/2-1);
    ny1=ceil(-n_new/2);ny2=ceil(n_new/2-1);
    NX=(nx1:1:nx2)';NY=(ny1:1:ny2)';% Column vector
   
    [KX0,KY0]=meshgrid(dky*NY,dkx*NX);
    %arbritary value of eta and xi from  Gauss-FFT
    eta=0.9306; xi=0.9306;
    KX=KX0+xi*dkx;
    KY=KY0+eta*dky;
    k=sqrt(KX.^2+KY.^2); 
    z2=data;
    %all delta fractions
    d_frac=0.05:0.05:0.75;
    for cnt=1:length(d_frac)
        
        %delta2=(max(z2(:))+min(z2(:)))*d_frac(cnt);
        delta2=max(z2(:))*d_frac(cnt);
        ss2=0;
        for n=1:30
            n=n-1;
            vv2=((-k).^n).*((z2-delta2).^(n+1)-(-delta2).^(n+1));
            ss2=ss2+vv2./(factorial(n+1));
            
        end
        val1=exp(-k.*delta2).*ss2;
        val2=(exp(-k.*z2)-1)./(-k);
        err(cnt)=norm(val2-val1);
    end
    
    [mn,id]=min(err);
    %delta=(max(z2(:))+min(z2(:)))*d_frac(id);
    delta=max(z2(:))*d_frac(id);
end
