function physioRDecoGUI()
    % Create the main figure window
    fig = uifigure('Name', 'Physio and R-Deco File Combiner', ...
                   'Position', [100 100 400 300], ...
                   'Resize', 'off');

    % Create labels and buttons for file selection
    uilabel(fig, 'Position', [20 250 150 20], 'Text', 'Select Physio File:');
    physioButton = uibutton(fig, 'push', ...
                           'Position', [180 250 200 22], ...
                           'Text', 'Browse...', ...
                           'ButtonPushedFcn', @(btn, event) selectPhysioFile(btn, fig));

    uilabel(fig, 'Position', [20 200 150 20], 'Text', 'Select R-Deco File:');
    rDecoButton = uibutton(fig, 'push', ...
                          'Position', [180 200 200 22], ...
                          'Text', 'Browse...', ...
                          'ButtonPushedFcn', @(btn, event) selectRDecoFile(btn, fig));

    uilabel(fig, 'Position', [20 150 150 20], 'Text', 'Select Output File:');
    outputButton = uibutton(fig, 'push', ...
                           'Position', [180 150 200 22], ...
                           'Text', 'Browse...', ...
                           'ButtonPushedFcn', @(btn, event) selectOutputFile(btn, fig));

    % Create Process button
    processButton = uibutton(fig, 'push', ...
                            'Position', [150 50 100 22], ...
                            'Text', 'Process', ...
                            'ButtonPushedFcn', @(btn, event) processFiles(btn, fig));

    % Store file paths in the figure's UserData
    fig.UserData = struct('physioPath', '', 'rDecoPath', '', 'outputPath', '');

    % Function to select physio file
    function selectPhysioFile(~, ~)
        [file, path] = uigetfile('*.mat', 'Select Physio File');
        if isequal(file, 0)
            return; % User canceled
        end
        fullPath = fullfile(path, file);
        fig.UserData.physioPath = fullPath;
        physioButton.Text = ['Physio: ' file];
    end

    % Function to select R-Deco file
    function selectRDecoFile(~, ~)
        [file, path] = uigetfile('*.mat', 'Select R-Deco File');
        if isequal(file, 0)
            return; % User canceled
        end
        fullPath = fullfile(path, file);
        fig.UserData.rDecoPath = fullPath;
        rDecoButton.Text = ['R-Deco: ' file];
    end

    % Function to select output file
    function selectOutputFile(~, ~)
        [file, path] = uiputfile('*.mat', 'Save Annotated Physio File As');
        if isequal(file, 0)
            return; % User canceled
        end
        fullPath = fullfile(path, file);
        fig.UserData.outputPath = fullPath;
        outputButton.Text = ['Output: ' file];
    end

    % Function to process files
    function processFiles(~, ~)
        physioPath = fig.UserData.physioPath;
        rDecoPath = fig.UserData.rDecoPath;
        outputPath = fig.UserData.outputPath;

        % Check if all paths are provided
        if isempty(physioPath) || isempty(rDecoPath) || isempty(outputPath)
            uialert(fig, 'Please select all files before processing.', 'Error');
            return;
        end

        % Check if files exist
        if ~exist(physioPath, 'file') || ~exist(rDecoPath, 'file')
            uialert(fig, 'One or more input files do not exist. Please check the paths.', 'Error');
            return;
        end

        try
            % Load the preprocessed physio data into a structure to avoid workspace conflicts
            physioData = load(physioPath);
            rDecoData = load(rDecoPath);

            % Extract variables from physioData (ensure they exist)
            requiredVars = {'RESP', 'RPIEZO', 'MRTRIG', 'STIMTRIG', 'name', 'PIEZOF', 'PIEZOD'};
            for var = requiredVars
                if ~isfield(physioData, var{1})
                    error(['Missing required variable in physio file: ' var{1}]);
                end
            end

            % Assign variables locally
            RESP = physioData.RESP;
            RPIEZO = physioData.RPIEZO;
            MRTRIG = physioData.MRTRIG;
            STIMTRIG = physioData.STIMTRIG;
            name = physioData.name;
            PIEZOF = physioData.PIEZOF;
            PIEZOD = physioData.PIEZOD;

            % Extract R-peak data from R-Deco
            if ~isfield(rDecoData.data, 'R_loc')
                error('R-peak data (R_loc) not found in R-Deco file.');
            end
            t = rDecoData.data.R_loc; t = t{1,1};
            r = ones(1, length(t));
            for j = 1:length(t)
                [Y, M, D, H, MN, S] = datevec(t(j));
                r(j) = H*3600 + MN*60 + S;
            end

            % Create a new physio struct to store all data
            physio = struct();
            physio.RESP = RESP;
            physio.RPIEZO = RPIEZO;
            physio.MRTRIG = MRTRIG;
            physio.STIMTRIG = STIMTRIG;
            physio.PIEZOF = PIEZOF;
            physio.PIEZOD = PIEZOD;
            physio.name = name;
            physio.piezoout = r;  % Add R-peak annotations (no task-specific field)

            % Save the combined data to the specified output file
            save(outputPath, 'physio');
            % Create a simple green checkmark icon (16x16 pixels for msgbox)
            greenIcon = zeros(16, 16, 3, 'uint8'); % RGB image
            greenIcon(:,:,2) = 255; % Set green channel to max (pure green)
            greenIcon(6:11, 3:8, :) = 0; % Black checkmark-like shape (simplified)
            greenIcon(6:11, 9:14, :) = 0; % Second part of checkmark

            % Display the message box
            msgbox(['Awesome! Your annotated physio data is saved to: ' outputPath], 'Success', 'custom', greenIcon);
            catch ME
            uialert(fig, ['Error processing files: ' ME.message], 'Error');
        end
    end
end
