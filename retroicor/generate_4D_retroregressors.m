function generate_4D_retroregressors(REGRESSORS,boldfile)
    %create 4D volumes for each of the regressors to pass into 3dTproject
    %inputs: REGRESSORS matrix generated during SMS retroicor, full path to bold file
    
    %load bold file 
    info = niftiinfo(boldfile);
    tmp = strsplit(boldfile,'.');
    tmp2 = tmp{1};
    
    xdim = info.ImageSize(1);
    ydim = info.ImageSize(2);
    nslices = info.ImageSize(3);
    maxvol = info.ImageSize(4);
    reg_info = info;
    reg_info.Datatype = 'double';
    
    nreg = size(REGRESSORS,2);
    ALL_REG = [];
    
    for r = 1:nreg
        reg = REGRESSORS(:,r,:);
        REG = zeros(xdim,ydim,nslices,maxvol);
        for t = 1:maxvol
            for s = 1:nslices
                REG(:,:,s,t) = repmat(reg(t,1,s),xdim,ydim);       
            end
        end
        ALL_REG.(strcat('REG',num2str(r))) = REG;
        niftiwrite(REG,strcat(tmp2,'_retreg-',num2str(r),'.nii.gz'),reg_info) %save regressor as nifti volume
        %need to save as double precision, not int...takes up a lot of space
    end
