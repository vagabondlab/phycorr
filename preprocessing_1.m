%% Main script to filter heartbeat signal and take first derivative
% Updated by Mario Murakami - March 2025
% This script is used for parsing physio .MAT files and then
% filtering/differentiating the pulse signal.
% Use this script when only 1 piezo was used
%% Define source and output folder
sourceFolder = '/autofs/cluster/vagabond/USERS/MARIO/Projects/7T/sourcedata/derivatives/retroicor/physiodata_raw_from_roberta';
outputFolder = '/autofs/cluster/vagabond/USERS/MARIO/Projects/7T/sourcedata/derivatives/retroicor/physio_data_preprocessed';

% Ensure output folder exists
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Get list of .mat files in source folder
matFiles = dir(fullfile(sourceFolder, '*.mat'));

%% Iterate through physio and filter
for i = 1:length(matFiles)
    %% Load the physio file
    filename = matFiles(i).name;
    filepath = fullfile(sourceFolder, filename);
    load(filepath);
    
    % Create identifier with task name
    s = regexprep(filename, '[^a-zA-Z0-9]', '_'); % Replace non-alphanumeric characters with underscores
    s = strrep(s, '.', '');  % Remove any periods
    
    % Ensure name starts with a letter (if it starts with a number)
    if isstrprop(s(1), 'digit')
        s = ['x' s];  % Prepend 'x' to make it a valid name
    end
    
    % Define separate variables for each channel
    RESP = data(datastart(1):dataend(1));
    RPIEZO = data(datastart(2):dataend(2));
    STIMTRIG = data(datastart(3):dataend(3));
    MRTRIG = data(datastart(4):dataend(4));
    name = s;
    
    % Calculate physio sampling rate
    fs = 1000;
    
    % Filter and differentiate piezo (when only 1). No LP channel used
    [PIEZOF, PIEZOD] = filter_hb(RPIEZO, fs);
    
    % Save filtered physio to new directory with unique filename
    outputFilename = [s '_physio.mat'];
    outputFilePath = fullfile(outputFolder, outputFilename);
    save(outputFilePath, 'RESP', 'RPIEZO', 'MRTRIG', 'STIMTRIG', 'name', 'PIEZOF', 'PIEZOD');
end

function [hbf, hbd] = filter_hb(hb, fs)
    %% Filter heartbeat signal then take first derivative
    % ---------------------------
    % INPUTS:
    % ---------------------------
    % hb (array): this is the raw heartbeat signal (ECG or pulse)
    % fs (scalar): sampling frequency in Hz
    % ---------------------------
    % OUTPUTS:
    % ---------------------------
    % hbf (array): filtered heartbeat signal
    % hbd (array): filtered/differentiated heartbeat signal
    
    % Define filter parameters
    f = [0 5 15 100] / (fs/2); % Adjusted frequency bands for stability
    a = [1 1 0 0];             % Passband: 0-5 Hz, Stopband: 15-100 Hz
    filterOrder = round(fs/2); % Increase filter order for smoother transition
    
    % Design the filter using firpm
    try
        h = firpm(filterOrder, f, a);
    catch ME
        warning('firpm failed to converge: %s. Switching to fir1 fallback.', ME.message);
        h = fir1(filterOrder, [5 15]/(fs/2), 'bandpass'); % Fallback to simpler bandpass filter
    end
    
    % Check for finite coefficients
    if ~all(isfinite(h))
        warning('Filter coefficients contain NaN or Inf. Using fallback filter.');
        h = fir1(filterOrder, [5 15]/(fs/2), 'bandpass'); % Fallback filter
    end
    
    % Apply filter
    hbf = filtfilt(h, 1, hb);
    
    % Take first derivative
    hbd = diff(hbf);
end
