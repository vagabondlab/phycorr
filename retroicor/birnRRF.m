function [RVT_RRF, TR_time] = birnRRF(TR,tot_time,RVT)

% TR and total acquisition time from fMRI data
% TR=3; tot_time=414*3; % in second
TR_time=[0:TR:tot_time-TR]';

% convolution function for respiratory signal (reference of Birn's paper)
RRF = @(t) 0.6*t^2.1*exp(-t/1.6) - 0.0023*t^3.54*exp(-t/4.25); 
Rrvt=[]; % response function
for nt=0:TR:40 % 0sec~28sec, interval: TR
    Rrvt=[Rrvt; RRF(nt)];
end

[p1, s1] = polyfit(TR_time,RVT,1);
lin_RVT = polyval(p1,[1:tot_time/TR]');  % detrending
RVTd = RVT-lin_RVT;

ConvRVT=conv(RVTd,Rrvt); % after convolution for RVT   
RVT_RRF=ConvRVT(1:tot_time/TR);