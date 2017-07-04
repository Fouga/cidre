function [cidre_correction,model] = check_cidreModel(coords,userConfig)

avDir = [userConfig.subdir.rawDataDir,filesep,userConfig.subdir.averageDir];
layer=coords(2); %optical section
chan=coords(5);
optical_section = 0;
% find mat files
fname = [avDir sprintf('/cidre_chanel%i_optical_section_%i.mat',chan,optical_section)];
% how cidre model structure should look like
model_def.method    = 'CIDRE';
model_def.v         = [];
model_def.z         = [];
model_def.v_small   = [];
model_def.z_small   = [];

if ~isempty(exist(fname,'file'))
    load(fname);
    names_def = fieldnames(model_def);
    names = fieldnames(model);
    tf = strcmp(names_def,names);
    if ~isempty(find(tf==0))
        cidre_correction = 0;
        disp('the fields in the cidre model are not correct.')
    else 
        cidre_correction = 1;
        
    end
else 
    cidre_correction = 0;
end