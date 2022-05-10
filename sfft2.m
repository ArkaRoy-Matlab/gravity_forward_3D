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

function X=sfft2(F,xi,eta,alpha,beta)
% 2d inverse shift fft of a spectrum F
% both X and F stored in natural orders
% xi,eta: frequency domain shift parameter of F along two dimensions
% alpha,beta: space domain shift parameter of X along two dimensions
%
F_Y=sfft_Y(F,eta,beta);
X=sfft_X(F_Y,xi,alpha);
