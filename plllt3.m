% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab code for plotting of rmse and time for different models for
%%different grid expansion ratio of quadrature based standard FFT method 
clear all
close all
%importing time data for different models for different grid expansion ratio
time=importdata(fullfile('.', 'output','sensitivity_time.txt'));
rmse=importdata(fullfile('.', 'output','sensitivity_rmse.txt'));

%figure for computation time for different grid expansion ratio
L=1:15;
figure(1)
hold on
for i=1:4
    t=time(:,i); 
    plot(L,t,'linewidth',2)
end
legend('Fixed density Model','Exponential density Model','Polynomial density Model','Complex density Model','location','best')
xlabel('Grid expansion ratio')
ylabel('Computational time (Sec)')
box on
%title('Plot for computational time due to standard FFT with grid expanded forward models vs. different grid expansion ratios.')

figure(2)
for i=1:4
    rm=rmse(:,i); 
    semilogy(L,rm,'linewidth',2)
    hold on
end
legend('Fixed density Model','Exponential density Model','Polynomial density Model','Complex density Model','location','best')
xlabel('Grid expansion ratio')
ylabel('log(RMS error)')
box on
%title('Plot for log(RMS error) between true gravity anomalies and gravity anomalies due to standard FFT with grid expanded forward models vs. different grid expansion ratios.')



   

