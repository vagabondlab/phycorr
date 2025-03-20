# retroicor

run in ordem, use Matlab
preprocessing_1.m, generate_1D_main_3v2.m, consolidation_2.m

The other script (generate_1D_fun_1.m) is a tool for the second script

## How this pipeline works?
First run the preprocessing_1.m using the raw physio data, then use R-DECO to extract the physio-peaks, 
Second run the consolidation_2.m to merge the files, in the GUI use the preprocess physio, the peaks extract from R-DECO and create a .mat file to output
Third run the generate_1D_main_3v2.m remember to addpath the generate_1D_fun.m, just select the folder that contains the merged files and the output file
