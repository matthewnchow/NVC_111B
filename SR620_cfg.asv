function [cmds, header] = SR620_cfg(cfg_filename, cfg_dir)
% Read configuration file and prepare list of commands
    switch nargin
        case 0
            cfg_filename = 'default_count.cfg'; %Put the name of the config file here
            cfg_dir = pwd;
        case 1
            cfg_dir = pwd;
    end
    
    cfg_filename = strcat(cfg_dir, '\SR620_configs\', cfg_filename);
    fid = fopen(cfg_filename);
    
    tline = fgetl(fid); %Header line
    i = 1;
    cmds = [];
    indices = [1];
    header = 'Output for: ';
    firstline = true; 
    while ischar(tline)
        if firstline
            header = strcat(header, tline, '\r\n');
            firstline = false;
        else
            cmds = [cmds tline];
            i = i + length(tline);
            indices = [indices, i];
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
end 