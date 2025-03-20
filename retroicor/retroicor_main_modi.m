function [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = ...
    retroicor_main_modi(image_matrix,slice_timing,TR,cardiac_trig_times,resp,OPTIONS)
%function [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = ...
%    retroicor_main(image_matrix,slice_timing,TR,cardiac_trig_times,resp,OPTIONS);
% 
% *  catie chang, 12/13/11
% *  adapted for more flexibile slice timing inputs  
%
% ---------------------------
% INPUTS:
% ---------------------------
% * image_matrix: 4D matrix of fmri data
% * slice_timing:  vector indicating slice acquisition times (within 1 TR), in seconds
%                  (the exact timing for physio phase calculations)
%                  (e.g., ([1 2 3 4 5 ... 30]-0.5)/30*TR for 30 slices)
% * TR: in seconds
% * cardiac_trig_times: vector of cardiac trigger times, in seconds.
% * resp: structure with field
%           * resp.trig: vector of respiratory trigger times, in sec.
%           -- OR instead with the following 2 fields -- 
%           * resp.wave: respiration amplitude signal   *AND*
%           * resp.dt:  sampling interval between the points in respiration
%            amplitude signal (in seconds, e.g. resp.dt=0.02 for 50 Hz sampling)
% * OPTIONS: optional structure with fields
%           * OPTIONS.slice_delta: constant temporal offset (in sec) to add
%           to all slice acquisition times  [default = 0]
%           * OPTIONS.doCorr:  =1 to correct the data, =0 if you
%           only want the regressors and phase outputs  [default =
%           1]
%           * OPTIONS.verbose: =1 to display progress messages
%           [default = 1]
%
%  (** set cardiac_trig_times = [] to do respiration only)
%  (** set resp = [] to do cardiac only)
%
%
% ---------------------------
% OUTPUTS:
% ---------------------------
% * image_matrix_corrected: 4D matrix of fmri data with "retroicor"
% terms regressed out
% * PHASES: cell array of cardiac & respiration phases for each
%      slice. PHASES{i}(:,1) contains the cardiac phase for
%      slice "i", and PHASES{i}(:,2) contains the resp phases for
%      slice "i".
% * REGRESSORS: the "retroicor" regressors for each
%      slice. dimensions = #timepoints x #regressors x # slices. The
%      matrix of regressors for slice "i" are REGRESSORS(:,:,i). 
% * OTHER: possibly useful stuff, see code
%
%


% defaults
delta = 0;
doCorr = 1;
verbose = 1;
CARD = 1;
RESP = 1;

% check options
if ~isempty(OPTIONS)
  if isfield(OPTIONS,'slice_delta')
    delta = OPTIONS.slice_delta;
  end
  if isfield(OPTIONS,'doCorr')
    doCorr = OPTIONS.doCorr;
  end
  if isfield(OPTIONS,'verbose')
    verbose = OPTIONS.verbose;
  end
end

if isempty(cardiac_trig_times)
  CARD = 0;
end
if isempty(resp)
  RESP = 0;
end
if (CARD+RESP==0)
  error('need resp and/or cardiac input')
end

% get image parameters
fmri_dims = size(image_matrix);
nslc = fmri_dims(3);
nframes = fmri_dims(4);
npix_x = fmri_dims(1);
npix_y = fmri_dims(2);

% cardiac input
if (CARD)
  etrig = cardiac_trig_times;
end

% respiration input
if (RESP)
  if isfield(resp,'trig')
    RESP_TRIGGERS = 1;
    rtrig = resp.trig;
  elseif isfield(resp,'wave')
    if ~isfield(resp,'dt')
      error('please specify resp.dt');
    end
    RESP_TRIGGERS = 0;
  else
    error('incorrect resp format');
  end
  
  if ~RESP_TRIGGERS
    % shift 
    respwave = resp.wave-min(resp.wave);
    % bin respiration signal into 100 values
    [Hb,bins] = hist(respwave,100);
    % calculate derivative
    % first, filter respiratory signal - just in case
    f_cutoff = 1; % max allowable freq
    fs = 1/resp.dt;
    wn = f_cutoff/(fs/2);
    ntaps = 20;
    b = fir1(ntaps,wn);
    respfilt = filtfilt(b,1,respwave);
    drdt = diff(respfilt); %derivative of filtered resp signal
  end
end

% find cardiac and respiratory phase vectors 
% (not yet accounting for slice acquisition order - that's next)

PHASES = {};
for jj=1:nslc
  % times at which ith slice was acquired (midpoint):
  slice_times = slice_timing(jj):TR:TR*nframes;
  slice_times = slice_times(1:nframes);
  % incorporate potential shift
  slice_times = slice_times + delta;
  
  phases_thisSlice = [];
  for ii=1:length(slice_times)
    
    % cardiac: find the closest R peaks on either side of the slice acquisition
        %phase is 2pi * (time since previous peak / RR interval length)
    if (CARD)
      prev_trigs = find(etrig<slice_times(ii));
      if isempty(prev_trigs)
        t1 = 0;
      else
        t1 = etrig(prev_trigs(end));
      end
      next_trigs = find(etrig>slice_times(ii));
      if isempty(next_trigs)
        t2 = nframes*TR;
      else
        t2 = etrig(next_trigs(1));
      end
      phi_cardiac = 2*pi*(slice_times(ii) - t1)/(t2-t1); %calculate phase
    else
      phi_cardiac = [];
    end
    
    if (RESP)
      if (RESP_TRIGGERS)
        % respiration: method based on triggers - same principle as cardiac
        % phase --> 2pi * (time into current cycle / resp cycle length)
        prev_trigs = find(rtrig<slice_times(ii));
        if isempty(prev_trigs)
          t1 = 0;
        else
          t1 = rtrig(prev_trigs(end));
        end
        next_trigs = find(rtrig>slice_times(ii));
        if isempty(next_trigs)
          t2 = nframes*TR;
        else
          t2 = rtrig(next_trigs(1));
        end
        phi_resp = (slice_times(ii) - t1)*2*pi/(t2-t1);
      else
        % respiration: method based on amplitude histogram
        tslice = slice_times(ii);
        iphys = max(1,round(tslice/resp.dt)); % closest index in resp waveform
        iphys = min(iphys,length(drdt)); %correct behavior at end of signal
        amp = respwave(iphys);
        dbins = abs(amp-bins);
        [blah,thisBin] = min(dbins);  %find closest resp histogram bin
        numer = sum(Hb(1:thisBin));
        phi_resp = pi*sign(drdt(iphys))*(numer/sum(Hb));
        %resp phases goes from -pi to pi because the histogram spans from
        %trough to peak, but can be inhalation or exhalation
        %taking the sum from 1:thisBin accounts for the the fact that it's
        %not a constant speed transition from peak to trough during
        %inhalation and exhalation
      end
    else
      phi_resp = [];
    end
    
    % store
    phases_thisSlice(ii,:) = [phi_cardiac, phi_resp];
  end
  
  
  PHASES{jj} = phases_thisSlice;
end

% generate regressors
REGRESSORS = [];
for jj=1:nslc
  
  if (CARD & RESP)
    phi_c = PHASES{jj}(:,1);
    phi_r = PHASES{jj}(:,2);
  elseif CARD
    phi_c = PHASES{jj}(:,1);
  elseif RESP
    phi_r = PHASES{jj}(:,1);
  end
  
  covs = [];
  % Fourier expansion of cardiac phase
  if (CARD)
    c1_c = cos(phi_c);
    s1_c = sin(phi_c);
    c2_c = cos(2*phi_c);
    s2_c = sin(2*phi_c);
    covs = [c1_c, s1_c, c2_c, s2_c];
  end
  % Fourier expansion of respiratory phase
  if (RESP)
    c1_r = cos(phi_r);
    s1_r = sin(phi_r);
    c2_r = cos(2*phi_r);
    s2_r = sin(2*phi_r);
    covs = [covs, c1_r, s1_r, c2_r, s2_r];
  end
  REGRESSORS(:,:,jj) = covs;
  
end


% correct the image data: slice-wise
PCT_VAR_REDUCED = zeros(npix_x,npix_y,nslc);
if (doCorr)
  image_matrix_corrected = zeros(size(image_matrix));
  for jj=1:nslc
    if verbose; fprintf('%d ... ',jj); end
    slice_data = squeeze(image_matrix(:,:,jj,:));
    Y_slice = (reshape(slice_data,npix_x*npix_y,nframes))'; %ntime x nvox
    Y_slice = cast(Y_slice,'double'); %convert to double precision for linear algebra operations
    t = [1:nframes]';
    XX = [t, t.^2, REGRESSORS(:,:,jj)];
    XX = [ones(size(XX,1),1), zscore(XX)]; %add columns of ones and normalize regressors
    if verbose; fprintf('regressing ... '); end
    Betas = pinv(XX)*Y_slice; %generate coefficients 
    Y_slice_corr = Y_slice - XX(:,4:end)*Betas(4:end,:); % take residuals but keep trends
    % calculate percent variance reduction (in voxel timeseries) using N-1 in var calculation 
    var_reduced = (var(Y_slice,0,1) - var(Y_slice_corr,0,1))./var(Y_slice,0,1);
    PCT_VAR_REDUCED(:,:,jj) = reshape(var_reduced',npix_x,npix_y); 
    % fill corrected volume
    V_slice_corr = Y_slice_corr';
    if verbose; fprintf('storing ... \n'); end
    for ii=1:nframes
      image_matrix_corrected(:,:,jj,ii) = reshape(V_slice_corr(:,ii),npix_x,npix_y); 
    end
  end
  fprintf('\n');
else
  % if you don't want to run the correction
  disp('NOT applying correction to data ... (returning regressors only).')
  image_matrix_corrected = image_matrix;
end

% return other possibly useful stuff
OTHER.PCT_VAR_REDUCED = PCT_VAR_REDUCED;
if (RESP)
  OTHER.drdt = drdt;
end
