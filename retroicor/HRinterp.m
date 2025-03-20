function [ hr, timevec ] = HRinterp( peaks, TR, maxvol )
%RRinterp returns the y values for given x values and plots RR differences
%and interpolation
%original sr / lambda = new sampling rate 
    % INPUTS:
        % peak locations in SECONDS!
        % TR and maxvol of scan (number of volumes)
    % OUTPUTS:
        % hr : heart rate linearly interpolated onto the TR locations
        % (timevec) 
    
hr_min=30; hr_max=180;

timevec = 0:TR:(maxvol*TR - TR);

rr = diff(peaks);
hr_inst = 60 ./ rr;
% check units
if mean(rr) >= 3
    error('peaks likely not in seconds!!')
end
hr_time = peaks(1:end-1) + rr / 2;  % midpoints between R waves

keep_inds = hr_inst <= hr_max & hr_inst >= hr_min;  % remove outliers 

hr_time = hr_time(keep_inds);
hr_inst = hr_inst(keep_inds);

hr = interp1(hr_time, hr_inst, timevec,'linear');

% replace nan values with mean HR
na_ind = isnan(hr);
hr(na_ind) = mean(hr,'omitnan');

%plot(hr_time, hr_inst)
%hold on
%plot(timevec, hr)
%title('interpolated RR peaks PPG')
end