%% PREPARE INPUTS FOR RETROICOR
% GUTBRAIN PROJECT
%written by Harrison Fisher 4/15/2020

%note on expected physio file format:
    % annoated physio data should be saved as a .mat file containing the
    % structure `Data` that has fields Filt_ECG (or PPG), ECG_Peaks (or PPG
    %_Peaks), and Trigger. Modify the code below if you have differently
    %named fields 

%INSTRUCTIONS
    %1) Set paths appropriately
    %2) Set Scan Parameters (vols, TR, sr)
    %3) Determine if correcting SMS squences or not
    %4) Input Subject identifier 
    %5) Run script:
        %a) Press any key after the first trigger plot shows up, if it
            %comes up again, it means that an incorrect number of triggers
            %were detected
        %b) On that second trigger plot, use the mouse to click on either
            %the start or the end of the triggers (whichever you are sure is
            %accurate)
        %c) if it is impossible to determine the exact start or stop, skip
            %this file and don't run retroicor on the corresponding bold data
    %6) pay attention to quality of the signals (ECG/PPG,
        %Resp, and trigger. Keep track of which are good enough to use during
        %the RETROICOR correction
physiodir = '/Users/hfisher/Dropbox (Partners HealthCare)/Pudding Study/data/';
annotdir = strcat(physiodir,'Annotated_MRI_Physio_2020.03./');
dest = strcat(physiodir,'GB_RETROICOR/');

vols = 288;
TR = 1.27;
sr = 500;
usePPG = 0;
fs = 40; %desired sampling frequency for the output files (default is 40Hz)

SMS=1; %set to 1 if generating files for SMS bold sequences, 0 if interleaved or ascending slice timing

%% Loop through subjects 

subj_list = {'FD_014'};
for j=1:numel(subj_list)
    subj = subj_list{j};
    list_raw = dir(strcat(annotdir,'*',subj,'*'));

    for i = 1:numel(list_raw)
        %i =1;
        load(strcat(annotdir,list_raw(i,1).name));
       
        %check to see if files have been created already    
        tmp = strsplit(list_raw(i).name,'.');
        tmp2 = strcat('RETRO_qrs_',tmp{1});
        if exist(strcat(dest,tmp2,'.1D'))==2 
            continue
        end
                
        % find triggers
        MRtrig = Data.Trigger; 
        [pks,locs] = findpeaks(MRtrig,'MINPEAKHEIGHT',0.5);
        locs_d = diff(locs);
        % remove extra trigger at the end of the signal if applicable
        if numel(locs) > vols && max(locs_d) > 1000
            [~, extra_trigger] = max(locs_d);
            locs(extra_trigger) = [];   pks(extra_trigger) = [];
        end
        %check trigger
        figure
        plot(Data.Trigger)
        hold on
        plot(locs,Data.Trigger(locs),'ro')
        hold off
        pause()
        close()
        
        %code to deal with improper trigger
        disp(numel(locs));
        start = locs(1);
        stop = locs(end);
        if numel(locs)~=vols
            plot(Data.Trigger)
            hold on
            plot(locs,Data.Trigger(locs),'ro')
            [x,~] = ginput(1);
            hold off
            close()

            if x>numel(Data.Trigger)/2
                %picked the end point
                [~, index] = min(abs(locs-x)); %find last trigger
                start = locs(index)-(vols-1)*TR*sr;
                stop = locs(index);
            else
                %picked starting point
                [~, index] = min(abs(locs-x));
                start = locs(index);
                stop = locs(index) + (vols-1)*TR*sr;
            end
            %generate "new" Trigger signal from start to stop
            locs = start:sr*TR:stop;
            %find actual trigger peaks for plotting 
            new_locs = correct_trigger_peaks(locs,Data.Trigger,vols,TR,sr);
            
            figure
            plot(Data.Trigger)
            hold on
            plot(new_locs,Data.Trigger(new_locs),'ro')
            hold off
            pause()
            close()
        end
        
        if usePPG==1
            RR = Data.PPG_Peaks(:,1);
            Signal = Data.PPG;
        else
            RR = Data.ECG_Peaks(:,1);
            Signal = Data.Filt_ECG;    
        end
        
        figure
        plot(Signal)
        hold on
        plot(RR,Signal(RR),'ro')
        plot([start,start],[min(Signal),max(Signal)],'g','linewidth',2) %start of scan
        plot([stop,stop],[min(Signal),max(Signal)],'r','linewidth',2) %end of scan
        hold off
        pause()
        close()

        %QRS 
        if SMS ==1
            %prepare cardiac trigger times in seconds from start of scan for SMS retroicor
            cut_RR = RR(RR >= start);
            RRretro = (cut_RR - start)/sr; 
        elseif SMS == 0
            %prepare binary signal of ones and zeros for afni 3dretroicor  
            QRS = zeros(1,numel(Signal));
            QRS(RR) = 1;     % qrs binary annotation
            QRS_cut = QRS(start:stop+TR*sr);       % cut pre- and post-imaging
        
            %need to smooth so that when you downsample you don't miss
            %peaks
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
          
        % respiration
        resp = Data.Resp(start:stop+TR*sr);
        temp = resp(1:round(sr/fs):end);
        RESPretro = medfilt1(temp,10);
        RESPretro(1) = temp(1);

        % plot
        figure
        subplot(2,1,1), plot((1:length(RESPretro))/fs,RESPretro);
        xlabel('Time (s)'); ylabel('Respiration')
        if SMS ==1
            subplot(2,1,2), plot(RRretro(2:end),diff(RRretro))
            xlabel('Time (s)'); ylabel('RR interval')
        elseif SMS ==0
            subplot(2,1,2), plot(1:length(RRretro),RRretro)
            xlabel('Time (s)'); ylabel('Binary QRS signal')
        end
        title(strcat(Name),'Interpreter', 'none');
        pause()
        close()

        % save
        fid = fopen(strcat(dest,'RETRO-qrs_',Name,'.1D'),'w');
        fprintf(fid,'%g\n',RRretro);
        fclose(fid);
        fid = fopen(strcat(dest,'RETRO-resp_',Name,'.1D'),'w');
        fprintf(fid,'%.4f\n',RESPretro);
        fclose(fid);

    end
end