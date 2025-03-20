function [RVHRcorr_reg] =  RVHR_corr_regressors(QRSretro_trig,respwave,TR,maxvol,fs,OPTIONS)
%RVHR_corr_SMS function to generate heart rate and respiratory volume per
%time regressors as well as  their first derivatives
%   slow fluctuations in physiological signals are not captured by
%   RETROICOR regressors that estimate noise components related to the
%   phases of cardiac and respiratory signals. Based on Chang and Glover
%   (2008, NeuroImage) and adapted by Harrison Fisher 

    % INPUTS
        % QRS_retro_trig: vector of cardiac trigger locations in SECONDS
        % respwave: respiratory signal (sampled at frequence fs)
        % TR: repetition time of bold data
        % maxvol: number of volumes in bold data
        % fs: sampling rate of respwave in Hz
        % OPTIONS: struct with fields savefile and filename
            % set savefile to 1 to save the regressors as a text file
            % specify name of text file in filename field
    % OUTPUTS
        % RVHRcorr_reg: convolved regressors sampled at each TR along with
        % the first derivatives
            
if isempty(OPTIONS)
    OPTIONS.savefile = 0;
end
            
colnames = {'RVconv','RVconv_d','HRconv','HRconv_d'};
%% HR processing and convolution
if isempty(QRSretro_trig)
    print('no QRS')
    hr_conv = [];
    hr_conv_d = [];
    colnames = {'RVconv','RVconv_d'};
else
    % interpolate & resample RR signal
    [hr, ~] = HRinterp(QRSretro_trig, TR, maxvol);
    
    % convolve HR with specific transfer function
    tot_time = TR*maxvol;
    hr_conv = zhangHRF(TR,tot_time,hr);
    if numel(hr_conv) > maxvol
        long = numel(hr_conv) - maxvol;
        hr_conv = hr_conv(long+1:end);
    end
    hr_conv_d = diff(hr_conv);
    hr_conv_d = [hr_conv_d(1); hr_conv_d];
end
%save?
%% Respiratory wave processing and convolution
if isempty(respwave)
    print('no respiration')
    rv_conv = [];
    rv_conv_d = [];
    colnames = {'HRconv','HRconv_d'};

else
    rv = RVTestimate(respwave,TR,maxvol,fs); % generate signal using st. dev. of resp wave over a 6s sliding window 
    rv_dm = rv - mean(rv);
    %convolve with respiratory response function
    tot_time = TR*maxvol;
    [rv_conv, ~] = birnRRF(TR,tot_time,rv_dm);
    if numel(rv_conv) > maxvol
        long = numel(rv_conv) - maxvol;
        rv_conv = rv_conv(long+1:end);
    end
    rv_conv_d = diff(rv_conv);
    rv_conv_d = [rv_conv_d(1); rv_conv_d];
end

RVHRcorr_reg = [rv_conv, rv_conv_d, hr_conv, hr_conv_d];

if isempty(RVHRcorr_reg)
    disp('no physio data')
else
    reg_table = array2table(RVHRcorr_reg);
    reg_table.Properties.VariableNames = colnames;
    if OPTIONS.savefile == 1
        writetable(reg_table,OPTIONS.filename)
    end
end
end
