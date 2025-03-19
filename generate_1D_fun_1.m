
function [qrs_1D_array, resp_1D_array] = generate_1D_fun_1(HB,R,resp,vols,...
        TR, sr, fs, SMS, MRtrig, fname)
%GENERATE_1D create .1D files to denoise fMRI data.
%   Use the respiration and heart beat signal (PPG or ECG) to generate 
%   1D files, which will be used to remove noise from respiration and blood flow.
%   The output file is saved to the current working directory. This
%   function also returns the 1D arrays.
% ---------------------------
% INPUTS:
% ---------------------------
% * HB (array): heart beat signal (PPG or filtered ECG)
% * R (array): peaks of the heartbeat signal in seconds.
% * resp (array): respiration signal
% * vols (double): number of volumes in scan  
% * TR (double): acqusition time (seconds)
% * sr (double): physio data sampling frequency
% * fs (double): Desired sampling frequency for output files (default = 40)
% * SMS (logical): 1 = generating files for SMS bold sequences 
%                  0 = interleaved or ascending slice timing
% * MRtrig (array): scanner trigger
% * fname (string): subject ID and other file identifiers
% ---------------------------
% OUTPUTS:
% ---------------------------
% * qrs_1D_array (array): QRS retroicor 1D array
% * resp_1D_array (array): respiration retroicor 1D array

%% Remove extra MR trigger if needed
% Take derivative of MR scanner trigger, removed diff
MRtrig_diff = diff(MRtrig);
% Find the peak of the derivative, which corresponds to peak of trigger
[~,trig_idx] = findpeaks(MRtrig_diff,'MINPEAKHEIGHT',2);

% Verify correct number of triggers
if length(trig_idx) ~= vols
    error('Number of MR scanner triggers ~= number of volumes.');
end

% Ensure the triggers are spaced apart appropriately (+/- 1 sample period)
dist = diff(trig_idx);
dmax = max(dist); t_dmax = dmax/sr;
dmin = min(dist); t_dmin = dmin/sr;
if t_dmax > TR + 1/sr
    error('Detected trigger peaks are too far apart.');
elseif t_dmin < TR - 1/sr
    error('Detected trigger peaks are too far apart.');
end

%% Plot trigger + detected trigger peaks. Assign start/stop indices
first_sample_idx = trig_idx(1);
last_sample_idx = (trig_idx(end)+TR*sr);

figure; hold on;
p(1) = plot(MRtrig);
p(2) = plot(trig_idx,MRtrig(trig_idx),'ro');
p(3) = xline(first_sample_idx,'g','linewidth',2); % start of sampling
p(4) = xline(last_sample_idx,'k','linewidth',2); % end of sampling
legend(p, {'Scanner Trigger','Scanner Trigger Annotated',...
    'Retroicor Sampling Start', 'Retroicor Sampling End'});
title({'Scanner Trigger with Annotated Peaks',...
    'Visually verify that the first/last trigger are correctly anotated.'});

%% Generate respiration retroicor array.
% filter resipration signal
resp = resp(first_sample_idx:last_sample_idx);
temp = resp(1:round(sr/fs):end);
RESPretro = medfilt1(temp,10);
RESPretro(1) = temp(1);

% save resp 1D
dest = strcat(pwd, '\');
fullfile = strcat(dest,'RETRO-resp_',fname,'.1D');
fid = fopen(char(fullfile),'w');
fprintf(fid,'%.4f\n',RESPretro);
fclose(fid);
resp_1D_array = RESPretro;

%% Convert the R-peaks into indices if not already
% Check if the array contains values to the right of the decimal.
b = mod(R(1),1);
% If so, then the R array is timestamps and not indices.
if b>0
    % Multiply by sampling rate to convert to indices
    R = uint64(R.*sr);
end

%% Generate QRS array
if ~isempty(R)
    %% Generate R retroicor array.
    if SMS ==1
        %prepare cardiac trigger times in seconds from start of scan for SMS retroicor
        cut_RR = R((first_sample_idx <= R)&(R <= last_sample_idx));
        cut_RR = double(cut_RR);
        % Subtract start index from R array. Divide by sampling rate.
        RRretro = (cut_RR - first_sample_idx)/sr; 
    elseif SMS == 0
        %prepare binary signal of ones and zeros for afni 3dretroicor  
        QRS = zeros(1,numel(HB));
        QRS(R) = 1;     % qrs binary annotation
        QRS_cut = QRS(first_sample_idx:last_sample_idx); % cut pre- and post-imaging

        %need to smooth so that when you downsample you don't miss peaks
        windowSize = 30; 
        b = (1/windowSize)*ones(1,windowSize); a = 1;
        QRS_sm = filter(b,a,QRS_cut);     % lp filtering with moving avg
        QRS_sm_rs = QRS_sm(1:round(sr/fs):end);      % 40Hz resampling
        [~,qrs] = findpeaks(QRS_sm_rs);
        RRretro = zeros(1,numel(QRS_sm_rs));
        RRretro(qrs) = 1;   
        RRretro = RRretro';
        if sum(QRS_cut) ~= sum(RRretro)
            disp('choose a different window size for smoothing') 
        end
        clear temp  
    end
    
    %% Plot sampled heartbeat with annotated peaks and trigger
    figure; hold on;
    p(1) = plot(HB);
    p(2) = plot(cut_RR,HB(cut_RR),'ro');
    p(3) = xline(first_sample_idx,'g','linewidth',2); % start of sampling
    p(4) = xline(last_sample_idx,'k','linewidth',2); % end of sampling
    legend(p, {'Heart Beat','Annotated Sampled Heart Beat Peaks',...
        'Retroicor Sampling Start', 'Retroicor Sampling End'});
    title({'Heart Beat, Annotated Beat, Retroicor Sampling Start/End',...
        'Visually verify that the heartbeat is correctly annotated.'})
    hold off;
    
    %% Plot respiration signal and R Interval
    figure;
    subplot(2,1,1), plot((1:length(RESPretro))/fs,RESPretro);
    xlabel('Time (s)'); ylabel('Respiration')
    if SMS ==1
        subplot(2,1,2), plot(RRretro(2:end),diff(RRretro))
        xlabel('Time (s)'); ylabel('R interval')
    elseif SMS ==0
        subplot(2,1,2), plot(1:length(RRretro),RRretro)
        xlabel('Time (s)'); ylabel('Binary QRS signal')
    end
    sgtitle('Respiration and R Interval')
    
    %% Save 1D file
    fullfile = strcat(dest,'RETRO-qrs_',fname,'.1D');
    fid = fopen(char(fullfile),'w');
    fprintf(fid,'%g\n',RRretro); fclose(fid);
    qrs_1D_array = RRretro;
else
    qrs_1D_array = [];
end

end
