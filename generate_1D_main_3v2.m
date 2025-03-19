function retroicorGUI()
    % Create the main figure window
    fig = uifigure('Name', 'RETROICOR 1D File Generator', ...
                   'Position', [100 100 400 350], ... % Increased height for new field
                   'Resize', 'off');

    % Create label and button for selecting the folder with annotated physio files
    uilabel(fig, 'Position', [20 300 150 20], 'Text', 'Select Physio Folder:');
    physioButton = uibutton(fig, 'push', ...
                           'Position', [180 300 200 22], ...
                           'Text', 'Browse...', ...
                           'ButtonPushedFcn', @(btn, event) selectPhysioFolder(btn, fig));

    % Create label and button for selecting the output folder
    uilabel(fig, 'Position', [20 250 150 20], 'Text', 'Select Output Folder:');
    outputButton = uibutton(fig, 'push', ...
                           'Position', [180 250 200 22], ...
                           'Text', 'Browse...', ...
                           'ButtonPushedFcn', @(btn, event) selectOutputFolder(btn, fig));

    % Create fields for user input (TR, sr)
    uilabel(fig, 'Position', [20 200 150 20], 'Text', 'TR (seconds):');
    trField = uieditfield(fig, 'numeric', ...
                          'Position', [180 200 100 22], ...
                          'Value', 1.19); % Default for MIGRAINE_PPG

    uilabel(fig, 'Position', [20 150 150 20], 'Text', 'Sampling Rate (Hz):');
    srField = uieditfield(fig, 'numeric', ...
                          'Position', [180 150 100 22], ...
                          'Value', 1000); % Updated to 1000 Hz as per your data

    % Removed SID field since it will be derived from filenames

    % Create Process button
    processButton = uibutton(fig, 'push', ...
                            'Position', [150 50 100 22], ...
                            'Text', 'Generate 1D Files', ...
                            'ButtonPushedFcn', @(btn, event) generate1DFiles(btn, fig));

    % Store data in the figure's UserData
    fig.UserData = struct('physioFolder', '', ...
                          'outputFolder', '', ...
                          'TR', 1.19, ...
                          'sr', 1000, ... % Updated to 1000 Hz
                          'fs', 40, ... % Constant for RETROICOR
                          'SMS', 1);   % Constant for SMS

    % Function to select physio folder
    function selectPhysioFolder(~, ~)
        folder = uigetdir('', 'Select Folder with Annotated Physio Files');
        if isequal(folder, 0)
            return; % User canceled
        end
        fig.UserData.physioFolder = folder;
        physioButton.Text = ['Folder: ' folder];
    end

    % Function to select output folder
    function selectOutputFolder(~, ~)
        folder = uigetdir('', 'Select Output Folder for 1D Files');
        if isequal(folder, 0)
            return; % User canceled
        end
        fig.UserData.outputFolder = folder;
        outputButton.Text = ['Output: ' folder];
    end

    % Function to generate 1D files
    function generate1DFiles(~, ~)
        % Retrieve values from UserData and fields
        physioFolder = fig.UserData.physioFolder;
        outputFolder = fig.UserData.outputFolder;
        TR = trField.Value;
        sr = srField.Value;
        fs = fig.UserData.fs; % Constant: 40 Hz
        SMS = fig.UserData.SMS; % Constant: 1

        % Check if folders are selected
        if isempty(physioFolder)
            uialert(fig, 'Please select a folder with annotated physio files.', 'Error');
            return;
        end
        if isempty(outputFolder)
            uialert(fig, 'Please select an output folder for the 1D files.', 'Error');
            return;
        end

        % Check if folders exist
        if ~exist(physioFolder, 'dir')
            uialert(fig, 'The selected physio folder does not exist.', 'Error');
            return;
        end
        if ~exist(outputFolder, 'dir')
            uialert(fig, 'The selected output folder does not exist.', 'Error');
            return;
        end

        % Get list of .mat files in the physio folder
        matFiles = dir(fullfile(physioFolder, '*.mat'));
        if isempty(matFiles)
            uialert(fig, 'No .mat files found in the selected folder.', 'Error');
            return;
        end

        try
            % Loop through each .mat file
            for i = 1:length(matFiles)
                % Full path to the current file
                physioPath = fullfile(physioFolder, matFiles(i).name);

                % Extract SID from the filename (remove extension)
                [~, sid, ~] = fileparts(matFiles(i).name);
                % Clean the SID (replace non-alphanumeric characters if needed)
                sid = regexprep(sid, '[^a-zA-Z0-9]', '_');
                if isstrprop(sid(1), 'digit')
                    sid = ['sub_' sid]; % Prepend to make it a valid identifier
                end

                % Load the physio struct from the .mat file
                loadedData = load(physioPath);
                if ~isfield(loadedData, 'physio')
                    warning('File %s does not contain a "physio" struct. Skipping...', matFiles(i).name);
                    continue;
                end
                physio = loadedData.physio;

                % Check required fields in physio struct
                requiredFields = {'RPIEZO', 'piezoout', 'RESP', 'MRTRIG'};
                for field = requiredFields
                    if ~isfield(physio, field{1})
                        warning('File %s is missing required field %s. Skipping...', matFiles(i).name, field{1});
                        continue;
                    end
                end

                % Debug: Check data types
                if ~isnumeric(physio.RPIEZO) || ~isnumeric(physio.piezoout) || ...
                   ~isnumeric(physio.RESP) || ~isnumeric(physio.MRTRIG)
                    warning('File %s has non-numeric fields. Skipping...', matFiles(i).name);
                    continue;
                end

                % Dynamically determine vols from MRTRIG
                MRtrig_diff = diff(physio.MRTRIG);
                [~, trig_idx] = findpeaks(MRtrig_diff, 'MINPEAKHEIGHT', 2);
                vols = length(trig_idx);
                if vols == 0
                    warning('File %s: No valid triggers detected in MRTRIG. Skipping...', matFiles(i).name);
                    continue;
                end

                % Generate 1D files with a base filename using the SID
                fnameBase = strcat('sub-', sid, '_run-001_bold');

                % Temporarily change the current directory to the output folder
                originalDir = pwd;
                cd(outputFolder);

                try
                    % Call generate_1D_fun_1
                    [~, ~] = generate_1D_fun_1(physio.RPIEZO, ... % HB: heart beat signal
                                              physio.piezoout, ... % R: R-peaks in seconds
                                              physio.RESP, ... % resp: respiration signal
                                              vols, ... % Number of volumes
                                              TR, ... % Acquisition time
                                              sr, ... % Sampling rate
                                              fs, ... % Desired sampling frequency
                                              SMS, ... % SMS flag
                                              physio.MRTRIG, ... % MR trigger
                                              fnameBase); % File identifier
                    fprintf('Successfully generated 1D files for %s\n', matFiles(i).name);
                catch ME
                    warning('Error processing file %s: %s. Skipping...', matFiles(i).name, ME.message);
                end

                % Restore the original directory
                cd(originalDir);
            end

            % Success message
            msgbox('RETROICOR 1D files generated successfully for all valid files!', 'Success');
        catch ME
            uialert(fig, ['Error generating 1D files: ' ME.message], 'Error');
        end
    end
end
