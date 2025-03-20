function retroicorGUI()
    % Create the main figure window
    fig = uifigure('Name', 'RETROICOR Processor', ...
                   'Position', [100 100 500 650], ...
                   'Resize', 'off');

    % Labels and buttons for file/folder selection
    uilabel(fig, 'Position', [20 600 150 20], 'Text', 'Select RETROICOR Folder:');
    retroFolderButton = uibutton(fig, 'push', ...
                                'Position', [180 600 200 22], ...
                                'Text', 'Browse...', ...
                                'ButtonPushedFcn', @(btn, event) selectRetroFolder(btn, fig));

    uilabel(fig, 'Position', [20 550 150 20], 'Text', 'Select BOLD File:');
    boldButton = uibutton(fig, 'push', ...
                         'Position', [180 550 200 22], ...
                         'Text', 'Browse...', ...
                         'ButtonPushedFcn', @(btn, event) selectBoldFile(btn, fig));

    uilabel(fig, 'Position', [20 500 150 20], 'Text', 'Select JSON File:');
    jsonButton = uibutton(fig, 'push', ...
                         'Position', [180 500 200 22], ...
                         'Text', 'Browse...', ...
                         'ButtonPushedFcn', @(btn, event) selectJsonFile(btn, fig));

    uilabel(fig, 'Position', [20 450 150 20], 'Text', 'Select QRS 1D File:');
    qrsButton = uibutton(fig, 'push', ...
                        'Position', [180 450 200 22], ...
                        'Text', 'Browse...', ...
                        'ButtonPushedFcn', @(btn, event) selectQRSFile(btn, fig));

    uilabel(fig, 'Position', [20 400 150 20], 'Text', 'Select RESP 1D File:');
    respButton = uibutton(fig, 'push', ...
                         'Position', [180 400 200 22], ...
                         'ButtonPushedFcn', @(btn, event) selectRESPFile(btn, fig));

    uilabel(fig, 'Position', [20 350 150 20], 'Text', 'Select Output Folder:');
    outputButton = uibutton(fig, 'push', ...
                           'Position', [180 350 200 22], ...
                           'Text', 'Browse...', ...
                           'ButtonPushedFcn', @(btn, event) selectOutputFolder(btn, fig));

    % Fields for TR and fs
    uilabel(fig, 'Position', [20 300 150 20], 'Text', 'TR (seconds):');
    trField = uieditfield(fig, 'numeric', ...
                          'Position', [180 300 100 22], ...
                          'Value', 1.19);

    uilabel(fig, 'Position', [20 250 150 20], 'Text', 'Sampling Frequency (Hz):');
    uilabel(fig, 'Position', [20 230 350 20], 'Text', 'Note: Set to match 1D file rate (default 40 Hz).');
    fsField = uieditfield(fig, 'numeric', ...
                          'Position', [180 250 100 22], ...
                          'Value', 40);

    % Optional features toggles
    uilabel(fig, 'Position', [20 200 150 20], 'Text', 'Optional Features:');
    saveHistCheck = uicheckbox(fig, 'Text', 'Save Histogram Plot', ...
                              'Position', [180 180 120 22], ...
                              'Value', 0);
    savePhaseCheck = uicheckbox(fig, 'Text', 'Save Phase Plot', ...
                               'Position', [180 150 120 22], ...
                               'Value', 0);
    saveVarNiiCheck = uicheckbox(fig, 'Text', 'Save Variance NIfTI', ...
                                'Position', [180 120 120 22], ...
                                'Value', 0);

    % Process button
    processButton = uibutton(fig, 'push', ...
                            'Position', [200 50 100 22], ...
                            'Text', 'Run RETROICOR', ...
                            'ButtonPushedFcn', @(btn, event) runRetroicor(btn, fig));

    % Store data in the figure's UserData
    fig.UserData = struct('retroFolder', '', ...
                          'boldPath', '', ...
                          'jsonPath', '', ...
                          'qrsPath', '', ...
                          'respPath', '', ...
                          'outputFolder', '', ...
                          'TR', 1.19, ...
                          'fs', 40, ...
                          'saveHist', 0, ...
                          'savePhase', 0, ...
                          'saveVarNii', 0);

    % Callback functions
    function selectRetroFolder(~, ~)
        folder = uigetdir('', 'Select RETROICOR Files Folder');
        if isequal(folder, 0)
            return;
        end
        fig.UserData.retroFolder = folder;
        retroFolderButton.Text = ['Folder: ' folder];
    end

    function selectBoldFile(~, ~)
        [file, path] = uigetfile({'*.nii;*.nii.gz', 'NIfTI Files (*.nii, *.nii.gz)'}, 'Select BOLD File', fig.UserData.retroFolder);
        if isequal(file, 0)
            return;
        end
        fullPath = fullfile(path, file);
        fig.UserData.boldPath = fullPath;
        boldButton.Text = ['BOLD: ' file];
    end

    function selectJsonFile(~, ~)
        [file, path] = uigetfile({'*.json', 'JSON Files (*.json)'}, 'Select JSON File', fig.UserData.retroFolder);
        if isequal(file, 0)
            return;
        end
        fullPath = fullfile(path, file);
        fig.UserData.jsonPath = fullPath;
        jsonButton.Text = ['JSON: ' file];
    end

    function selectQRSFile(~, ~)
        [file, path] = uigetfile({'*.1D', '1D Files (*.1D)'}, 'Select QRS 1D File', fig.UserData.retroFolder);
        if isequal(file, 0)
            return;
        end
        fullPath = fullfile(path, file);
        fig.UserData.qrsPath = fullPath;
        qrsButton.Text = ['QRS: ' file];
    end

    function selectRESPFile(~, ~)
        [file, path] = uigetfile({'*.1D', '1D Files (*.1D)'}, 'Select RESP 1D File', fig.UserData.retroFolder);
        if isequal(file, 0)
            return;
        end
        fullPath = fullfile(path, file);
        fig.UserData.respPath = fullPath;
        respButton.Text = ['RESP: ' file];
    end

    function selectOutputFolder(~, ~)
        folder = uigetdir('', 'Select Output Folder');
        if isequal(folder, 0)
            return;
        end
        fig.UserData.outputFolder = folder;
        outputButton.Text = ['Output: ' folder];
    end

    function runRetroicor(~, ~)
        % Retrieve values from UserData and fields
        retroFolder = fig.UserData.retroFolder;
        boldPath = fig.UserData.boldPath;
        jsonPath = fig.UserData.jsonPath;
        qrsPath = fig.UserData.qrsPath;
        respPath = fig.UserData.respPath;
        outputFolder = fig.UserData.outputFolder;
        TR = trField.Value;
        fs = fsField.Value;
        saveHist = saveHistCheck.Value;
        savePhase = savePhaseCheck.Value;
        saveVarNii = saveVarNiiCheck.Value;

        % Validation
        if isempty(retroFolder) || isempty(boldPath) || isempty(jsonPath) || isempty(outputFolder)
            uialert(fig, 'Please select RETROICOR folder, BOLD file, JSON file, and output folder.', 'Error');
            return;
        end
        if ~exist(boldPath, 'file')
            uialert(fig, 'The selected BOLD file does not exist.', 'Error');
            return;
        end
        if ~exist(jsonPath, 'file')
            uialert(fig, 'The selected JSON file does not exist.', 'Error');
            return;
        end
        if ~exist(retroFolder, 'dir')
            uialert(fig, 'The selected RETROICOR folder does not exist.', 'Error');
            return;
        end
        if ~exist(outputFolder, 'dir')
            uialert(fig, 'The selected output folder does not exist.', 'Error');
            return;
        end

        try
            % Test if the BOLD file is readable by niftiread
            try
                niftiread(boldPath);
            catch
                uialert(fig, 'The selected BOLD file is not a valid NIfTI file.', 'Error');
                return;
            end

            % Extract filename components
            [datafolder, funcnameWithExt, ~] = fileparts(boldPath);
            [funcname, ~] = fileparts(funcnameWithExt);
            outpath = fullfile(outputFolder, funcname);

            % Check if already processed
            outfunc = [outpath '_retro-corrected'];
            if exist([outfunc '.nii.gz'], 'file')
                uialert(fig, [funcname ' already corrected.'], 'Warning');
                return;
            end

            % Load physio files
            resp_struct = [];
            if ~isempty(respPath) && exist(respPath, 'file')
                RESPretro = readtable(respPath, 'FileType', 'text', 'ReadVariableNames', false);
                resp_struct.wave = RESPretro.Var1;
                resp_struct.dt = 1 / fs;
            else
                uialert(fig, 'No respiration file found or selected.', 'Warning');
            end

            QRSretro_trig = [];
            if ~isempty(qrsPath) && exist(qrsPath, 'file')
                QRStmp = readtable(qrsPath, 'FileType', 'text', 'ReadVariableNames', false);
                QRSretro = QRStmp.Var1;
                if max(QRSretro) == 1
                    QRSretro_trig = find(QRSretro == 1) / fs;
                elseif (1 / mean(diff(QRSretro))) > 20
                    QRSretro_trig = QRSretro / fs;
                else
                    QRSretro_trig = QRSretro;
                end
            else
                uialert(fig, 'No QRS file found or selected.', 'Warning');
            end

            % Load BOLD and JSON
            image_matrix = niftiread(boldPath);
            info = niftiinfo(boldPath);
            [xdim, ydim, nslices, maxvol] = size(image_matrix);
            if exist(jsonPath, 'file')
                js = jsondecode(fileread(jsonPath));
                ST = js.SliceTiming;
                if isfield(js, 'RepetitionTime')
                    TR = js.RepetitionTime; % Override TR if in JSON
                end
                js.CardiacPhysio = ~isempty(qrsPath) && exist(qrsPath, 'file');
                js.RespPhysio = ~isempty(respPath) && exist(respPath, 'file');
                jsonout = [outpath '.json'];
                fid = fopen(jsonout, 'w');
                fprintf(fid, jsonencode(js));
                fclose(fid);
            else
                uialert(fig, 'JSON file not found.', 'Error');
                return;
            end

            % Run RETROICOR
            disp(['Running RETROICOR for ' funcname]);
            OPTIONS.doCorr = 1;
            [image_matrix_corrected, PHASES, REGRESSORS, OTHER] = retroicor_main_modi(image_matrix, ST, TR, QRSretro_trig, resp_struct, OPTIONS);

            % Save corrected image
            image_matrix_corrected = int16(image_matrix_corrected);
            niftiwrite(image_matrix_corrected, [outfunc '.nii.gz'], info);

            % Save regressors and variance
            save([outpath '_retro-regressors.mat'], 'REGRESSORS');
            save([outpath '_retro-pctvar.mat'], '-struct', 'OTHER', 'PCT_VAR_REDUCED');

            % Generate RVHRcorr regressors
            if ~isempty(qrsPath) && exist(qrsPath, 'file') && ~isempty(respPath) && exist(respPath, 'file')
                RVHRoptions = [];
                RVHRoptions.savefile = 0; % We'll save HRconv separately
                RVHRcorr_reg = RVHR_corr_regressors(QRSretro_trig, resp_struct.wave, TR, maxvol, fs, RVHRoptions);
                if ~isempty(RVHRcorr_reg)
                    hr_conv = RVHRcorr_reg(:, 3); % Extract HRconv (column 3)
                    writematrix(hr_conv, [outpath '_RVHRcorr_HRconv_feat.txt'], 'Delimiter', ' ');
                    disp(['Saved HR-convolved regressor: ' outpath '_RVHRcorr_HRconv_feat.txt']);
                else
                    disp('No RVHRcorr regressors generated due to missing physio data.');
                end
            end

            % Save phases if both files exist
            if ~isempty(respPath) && ~isempty(qrsPath) && exist(respPath, 'file') && exist(qrsPath, 'file')
                card_phases = zeros(maxvol, nslices);
                resp_phases = zeros(maxvol, nslices);
                for sl = 1:nslices
                    phases = PHASES{sl};
                    card_phases(:, sl) = phases(:, 1);
                    resp_phases(:, sl) = phases(:, 2);
                end
                writetable(array2table(card_phases), [outpath '_retro-cardphases.txt']);
                writetable(array2table(resp_phases), [outpath '_retro-respphases.txt']);
            elseif ~isempty(respPath) && exist(respPath, 'file')
                resp_phases = zeros(maxvol, nslices);
                for sl = 1:nslices
                    phases = PHASES{sl};
                    resp_phases(:, sl) = phases(:, 1);
                end
                writetable(array2table(resp_phases), [outpath '_retro-respphases.txt']);
            elseif ~isempty(qrsPath) && exist(qrsPath, 'file')
                card_phases = zeros(maxvol, nslices);
                for sl = 1:nslices
                    phases = PHASES{sl};
                    card_phases(:, sl) = phases(:, 1);
                end
                writetable(array2table(card_phases), [outpath '_retro-cardphases.txt']);
            end

            % Optional outputs
            if saveVarNii
                var_info = info;
                var_info.PixelDimensions = info.PixelDimensions(1:3);
                var_info.ImageSize = info.ImageSize(1:3);
                pctvar = int16(round(OTHER.PCT_VAR_REDUCED * 10000));
                niftiwrite(pctvar, [outpath '_retro-PctVarReduced.nii.gz'], var_info);
            end

            if saveHist
                [~, EDGES] = histcounts(squeeze(OTHER.PCT_VAR_REDUCED(:,:,:)));
                adj_edges = EDGES(1):((EDGES(end)-EDGES(1))/50):EDGES(end);
                hist_mat = zeros(numel(adj_edges)-1, nslices);
                for sl = 1:nslices
                    Hb = histcounts(squeeze(OTHER.PCT_VAR_REDUCED(:,:,sl)), adj_edges);
                    hist_mat(:, sl) = Hb;
                end
                xx = adj_edges(1:end-1) + diff(adj_edges);
                S = figure('visible', 'off');
                surf(xx * 100, 1:nslices, hist_mat', 'FaceAlpha', 0.7);
                xlabel('% Variance Change'); ylabel('Slices'); zlabel('# voxels');
                colorbar;
                saveas(S, [outpath '_retro-varchgsurf.png']);
                close(S);
            end

            if savePhase
                if ~isempty(respPath) && ~isempty(qrsPath) && exist(respPath, 'file') && exist(qrsPath, 'file')
                    ph = figure('visible', 'off');
                    subplot(2,1,1);
                    heatmap(card_phases', 'Colormap', parula);
                    xlabel('Time (TRs)'); ylabel('Slices'); title('Cardiac Phase');
                    ax = gca; ax.XDisplayLabels = repmat(' ', maxvol, 1); ax.YDisplayLabels = repmat(' ', nslices, 1);
                    subplot(2,1,2);
                    heatmap(resp_phases' + pi, 'Colormap', parula);
                    xlabel('Time (TRs)'); ylabel('Slices'); title('Respiratory Phase');
                    ax = gca; ax.XDisplayLabels = repmat(' ', maxvol, 1); ax.YDisplayLabels = repmat(' ', nslices, 1);
                    saveas(ph, [outpath '_retro-phases.png']);
                    close(ph);
                elseif ~isempty(respPath) && exist(respPath, 'file')
                    ph = figure('visible', 'off');
                    heatmap(resp_phases' + pi, 'Colormap', parula);
                    xlabel('Time (TRs)'); ylabel('Slices'); title('Respiratory Phase');
                    ax = gca; ax.XDisplayLabels = repmat(' ', maxvol, 1); ax.YDisplayLabels = repmat(' ', nslices, 1);
                    saveas(ph, [outpath '_retro-phases.png']);
                    close(ph);
                elseif ~isempty(qrsPath) && exist(qrsPath, 'file')
                    ph = figure('visible', 'off');
                    heatmap(card_phases', 'Colormap', parula);
                    xlabel('Time (TRs)'); ylabel('Slices'); title('Cardiac Phase');
                    ax = gca; ax.XDisplayLabels = repmat(' ', maxvol, 1); ax.YDisplayLabels = repmat(' ', nslices, 1);
                    saveas(ph, [outpath '_retro-phases.png']);
                    close(ph);
                end
            end

            msgbox(['RETROICOR processing completed for ' funcname '!'], 'Success');
        catch ME
            uialert(fig, ['Error running RETROICOR: ' ME.message], 'Error');
        end
    end
end
