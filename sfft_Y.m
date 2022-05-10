% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% *** inspired by from the paper
%%"Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68"
% ****************************************************************
function F_Y=sfft_Y(F,eta,beta)
% 1d inverse shift fft along each column of F
% eta: frequency domain shift parameter
% beta: space domain shift parameter
%
[ny,nx]=size(F);
ny1=ceil(-ny/2);
ny2=ceil(ny/2-1);
NY=(ny1:1:ny2)';

wy=exp(-1i*2*pi/ny);
TR_y=wy.^(eta*(NY+beta));
TR_Y=repmat(TR_y,[1,nx]);
TL_y=wy.^(beta*NY);
TL_Y=repmat(TL_y,[1,nx]);
F1=ifftshift(F.*TR_Y);
X1=fft(F1);
F_Y=TL_Y.*fftshift(X1);