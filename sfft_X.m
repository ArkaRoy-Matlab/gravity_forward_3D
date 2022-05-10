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
function F_X=sfft_X(F,xi,alpha)
% 1d inverse shift fft along each row of F
% xi: frequency domain shift parameter
% alpha: space domain shift parameter
%
F=transpose(F);
F_Y=sfft_Y(F,xi,alpha);
F_X=transpose(F_Y);
