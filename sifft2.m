%Copyright (c) 2018, LeyuanWu
%Copyright (c) 2009, Greg von Winckel
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:

%* Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.

%* Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% *** Originated by:
% *** from the paper
%%"Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68"
% ****************************************************************

%%Gauss-FFT code to reproduce the numerical results in the paper
%%"Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68"
function X=sifft2(F,xi,eta,alpha,beta)
% 2d inverse shift fft of a spectrum F
% both X and F stored in natural orders
% xi,eta: frequency domain shift parameter of F along two dimensions
% alpha,beta: space domain shift parameter of X along two dimensions
%
F_Y=sifft_Y(F,eta,beta);
X=sifft_X(F_Y,xi,alpha);
