function [RVHR_HR] = zhangHRF(TR,tot_time,RR)

% TR and total acquisition time from fMRI data
% TR=3; tot_time=414*3; % in second
TR_time=[TR:TR:tot_time]';

% convolution function for cardiac signal (reference of Zhang's paper)
CRF = @(t) 0.6*t^2.7*exp(-t/1.6) - 16/sqrt(2*pi*9)*exp(-(t-12)^2/2/9); 
Rhr=[]; % response function
for nt=0:TR:28 % 0sec~28sec, interval: TR
    Rhr=[Rhr; CRF(nt)];
end

% HR = 60*ones(size(RR))./RR; % bpm
HR = RR;

if find(size(HR) > 1 == 1) == 2
    HR = HR';
end

[p1, s1] = polyfit(TR_time,HR,1);
lin_HR = polyval(p1,[0:numel(HR)-1]');  % detrending
HR = HR - lin_HR;

ConvHR=conv(HR,Rhr); % after convolution for HR   
RVHR_HR=ConvHR(1:tot_time/TR+1);