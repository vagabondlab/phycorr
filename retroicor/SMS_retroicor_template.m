%% RETROICOR FOR SMS BOLD DATA
%written by Harrison Fisher 4/16/2020 (during quarantine)

% note: script assumes that bold files end in bold.nii.gz and physio files having naming convention
%  RETRO-qrs/resp_"boldfilename".1D
    %e.g. for sub-GUTBRAINHC002_ses-MRI2A_task-rest_run-001_bold.nii.gz:
        %RETRO-qrs_sub-GUTBRAINHC002_ses-MRI2A_task-rest_run-001_bold.1D'
    %need to change code if your files are named differently 

%if QRS 1D file is peak locations, it should be in SECONDS. If it is in
%samples, the code will fix it, but safer to always make it in seconds. 

tempdir = '/autofs/space/shasta_002/users/SPARC_BRAINGUT/tmp/';  % set so that matlab doesn't run into space issues unzipping files in /tmp/

root = '/autofs/space/shasta_002/users/SPARC_BRAINGUT/';
bolddir = strcat(root,'GUTBRAIN_BIDS/');
physiodir = strcat(root,'physio/RETROICOR/'); %directory with RETROICOR input files 
rvhrcorrdir = strcat(root,'physio/RVHRcorr/'); %directory for RVHRcorr outputs
outdir = strcat(bolddir,'derivatives/retroicor/');  % directory for output files (should be "derivatives/retroicor" if in BIDS format
OPTIONS = [];
OPTIONS.doCorr = 1; % set to zero to skip applying the correction to the data and just generate the regressors

isBIDS = 1;
fs = 40; %sampling frequency of retroicor files
doRVHRcorr = 1 
TR = 1.27;

save_hist_plot = 0; %set to 1 to save 3d histogram surface of percent variance change
save_phase_plot = 0; %set to 1 to save heatmaps of physio signal phase timeseries
save_var_nii = 0; %set to 1 to 3D nii percent change in variance for voxel timeseries
    % note that values are the percentages are multiplied by 100 in order
    % to save as integers and reduce file size. 

subj_list = {'GUTBRAINHC'};
for j=1:numel(subj_list)
    subj = subj_list{j};
    list_raw = dir(strcat(bolddir,'*',subj,'*/ses-*/func/*bold.nii.gz'));

    for i = 1:numel(list_raw)

        funcname = list_raw(i).name;
        datafolder = strcat(list_raw(i).folder,'/');
        tmp = strsplit(funcname,'.');
        tmp2 = tmp{1};
        tmp_out = strsplit(datafolder,bolddir);
        inpath = strcat(datafolder,tmp2);
        mkdir(strcat(outdir,tmp_out{2}))
        outpath = strcat(outdir,tmp_out{2},tmp2);
              
        
        outfunc = strcat(outpath,'_retro-corrected.nii.gz'); %output file
        
        if exist(outfunc,'file') == 2
           disp(strcat(funcname,' already corrected'))
           continue
        end
        qrsfile = strcat(physiodir,'RETRO-qrs_',tmp2,'.1D');
        respfile =  strcat(physiodir,'RETRO-resp_',tmp2,'.1D');
        
        %skip to next iteration if no retroicor files
        if ~exist(respfile,'file') && ~exist(qrsfile,'file')
            disp(['NO PHYSIO DATA FOR ',tmp2])
            continue  
        end
        
        % Load Physio Files
        if exist(respfile,'file') == 2
            RESPretro = readtable(respfile,'FileType','text','ReadVariableNames',false);
            resp_struct = [];
            resp_struct.wave = RESPretro.Var1;
            resp_struct.dt = 1/fs;
        else
            resp_struct = [];
            disp(strcat(['NO RESPIRATION FILE FOUND FOR ',tmp2]))
        end

        if exist(qrsfile,'file') == 2
            QRStmp = readtable(qrsfile,'FileType','text','ReadVariableNames',false);
            QRSretro = QRStmp.Var1;
            if max(QRSretro) == 1
                %recover peak locations from binary qrs signal
                QRSretro_trig = find(QRSretro ==1)/fs;
            elseif (1/mean(diff(QRSretro))) > 20
                %units are samples not seconds...fix
                QRSretro_trig = QRSretro / fs;
            else
                QRSretro_trig = QRSretro;
            end   
        else
            QRSretro_trig = [];
            disp(strcat(['NO QRS FILE FOUND FOR ',tmp2]))
        end

        if isBIDS == 1
            %read slice timing from json file
            jsonname = strcat(inpath,'.json');
            js = jsondecode(fileread(jsonname));
            ST = js.SliceTiming;
            TR = js.RepetitionTime;
            if exist(qrsfile,'file') == 2
               js.CardiacPhysio=1;
            else
               js.CardiacPhysio=0;
            end
             if exist(respfile,'file') == 2
               js.RespPhysio=1;
            else
               js.RespPhysio=0;
            end
            % write physio info to json file
            jsonout = strcat(outpath,'.json');
            fid=fopen(jsonout,'w');
            encodedJSON = jsonencode(js); 
            fprintf(fid, encodedJSON); 
        else
           STfile = strcat(datafolder,'STfile.txt');
           st_tmp = readtable(STfile,'ReadVariableNames',false);
           ST = st_tmp.Var1;
           clear st_tmp
        end    
        
        image_matrix = niftiread(strcat(datafolder,funcname));
        info = niftiinfo(strcat(datafolder,funcname));
        [xdim,ydim,nslices,maxvol] = size(image_matrix);
        
        disp(strcat('running retroicor for ',funcname))
        
        if doRVHRcorr == 1
            %generate RVHRcorr Regressors 
            RVHRoptions=[];
            RVHRoptions.savefile=1;
            RVHRoptions.filename = strcat(rvhrcorrdir,tmp2,'_RVHRcorr.txt');
            RVHR_corr_regressors(QRSretro_trig,resp_struct.wave,TR,maxvol,fs,RVHRoptions);
        end
            
        %estimate phases, create regressors, and correct the image
        [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = retroicor_main_modi(image_matrix,ST,TR,QRSretro_trig,resp_struct,OPTIONS);
        
        if OPTIONS.doCorr == 1
            %save corrected image 
            image_matrix_corrected = int16(image_matrix_corrected);
            niftiwrite(image_matrix_corrected,outfunc,info)
        end

        %save matlab array of REGRESSORS
        save(strcat(outpath,'_retro-regressors.mat'),'REGRESSORS')
             
        %save matlab array of PCT variance change 
      	save(strcat(outpath,'_retro-pctvar.mat'),'-struct','OTHER','PCT_VAR_REDUCED')
      	
        %reformat PHASES into arrays
        if exist(respfile,'file') == 2 && exist(qrsfile,'file') ==2           
            card_phases = zeros(maxvol,nslices);
            resp_phases = zeros(maxvol,nslices);
            for sl = 1:nslices
               phases = PHASES{sl};
               card_phases(:,sl) = phases(:,1);
               resp_phases(:,sl) = phases(:,2);
            end
            %save phases in nvol x nslice table 
            writetable(array2table(card_phases),strcat(outpath,'_retro-cardphases.txt'))
            writetable(array2table(resp_phases),strcat(outpath,'_retro-respphases.txt'))         
        elseif exist(respfile,'file') == 2  % if RESP only
            resp_phases = zeros(maxvol,nslices);
            for sl = 1:nslices
               phases = PHASES{sl};        
               resp_phases(:,sl) = phases(:,1);
            end
            writetable(array2table(resp_phases),strcat(outpath,'_retro-respphases.txt'))         
        elseif exist(qrsfile,'file') == 2   % if QRS only
            card_phases = zeros(maxvol,nslices);
            for sl = 1:nslices
               phases = PHASES{sl};
               card_phases(:,sl) = phases(:,1);
            end
            writetable(array2table(card_phases),strcat(outpath,'_retro-cardphases.txt'))     
        end
        %% Extra Stuff
        if save_var_nii ==1
            var_info = info;
            var_info.PixelDimensions = info.PixelDimensions(1:3);
            var_info.ImageSize = info.ImageSize(1:3);
            pctvar = int16(round(OTHER.PCT_VAR_REDUCED*10000));
            niftiwrite(pctvar,strcat(outpath,'_retro-PctVarReduced.nii.gz'),var_info)
        end
        
        if save_hist_plot == 1
            %PCT_VAR_REDUCED is 3D image where each voxel is 
                % (var corrected - var orig) / var orig 

            [~,EDGES] = histcounts(squeeze(OTHER.PCT_VAR_REDUCED(:,:,:)));
            adj_edges = EDGES(1):((EDGES(end)-EDGES(1))/50):EDGES(end);
            hist_mat = zeros(numel(adj_edges)-1,nslices);
            for sl = 1:nslices
                Hb = histcounts(squeeze(OTHER.PCT_VAR_REDUCED(:,:,sl)),adj_edges);
                hist_mat(:,sl) = Hb;
            end
            xx = adj_edges(1:end-1) + diff(adj_edges);
            S = figure('visible', 'off');
            surf(xx*100,1:nslices,hist_mat','FaceAlpha',0.7);
            xlabel("% Variance Change"); ylabel("Slices");zlabel("# voxels")
            colorbar
            xh = get(gca,'XLabel'); % Handle of the x label
            set(xh, 'Units', 'Normalized')
            pos = get(xh, 'Position');
            set(xh, 'Position',pos.*[1,-0.5,1],'Rotation',15)

            yh = get(gca,'YLabel'); % Handle of the x label
            set(yh, 'Units', 'Normalized')
            pos = get(yh, 'Position');
            set(yh, 'Position',pos.*[1,-0.5,1],'Rotation',-30)
            
            saveas(S,strcat(outpath,'_retro-varchgsurf.png'))
            close(S)
            
        end
     
        if save_phase_plot == 1
            if exist(respfile,'file') == 2 && exist(qrsfile,'file') == 2       
                ph = figure('visible', 'off');
                subplot(2,1,1)
                heatmap(card_phases','Colormap',parula)
                xlabel("Time (TRs)");ylabel("Slices");title('Cardiac Phase')
                ax = gca;
                ax.XDisplayLabels = repmat([' '],maxvol,1);
                ax.YDisplayLabels = repmat([' '],nslices,1);
                subplot(2,1,2)
                heatmap(resp_phases' + pi,'Colormap',parula)
                xlabel("Time (TRs)"); ylabel("Slices");title('Respiratory Phase')
                ax = gca;
                ax.XDisplayLabels = repmat([' '],maxvol,1);
                ax.YDisplayLabels = repmat([' '],nslices,1);
                saveas(ph,strcat(outpath,'_retro-phases.png'))
                close(ph)
            elseif exist(respfile,'file') == 2
                ph = figure('visible', 'off');
                heatmap(resp_phases' + pi,'Colormap',parula)
                xlabel("Time (TRs)"); ylabel("Slices");title('Respiratory Phase')
                ax = gca;
                ax.XDisplayLabels = repmat([' '],maxvol,1);
                ax.YDisplayLabels = repmat([' '],nslices,1);
                saveas(ph,strcat(outpath,'_retro-phases.png'))
                close(ph)
            elseif exist(qrsfile,'file') == 2 
                ph = figure('visible', 'off');
                heatmap(card_phases','Colormap',parula)
                xlabel("Time (TRs)");ylabel("Slices");title('Cardiac Phase')
                ax = gca;
                ax.XDisplayLabels = repmat([' '],maxvol,1);
                ax.YDisplayLabels = repmat([' '],nslices,1);
                saveas(ph,strcat(outpath,'_retro-phases.png'))
                close(ph)
            end
        end

    end
    
end
