function SR620_Master()

    % General GPIB control script for SR620
    % Configures by input commands or a configuration file
    % Takes .cfg configuration file, which contains commands for SR620
    % Output of SR620 measurement written to a .txt file

    cmds = SR620_cfgfileread()

    % Find a GPIB object. Call it SR620
    sr620 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');

    % Create the GPIB object if it does not exist
    % otherwise use the object that was found.
    if isempty(sr620)
        sr620 = gpib('NI', 0, 16);
    else
        fclose(sr620);
        sr620 = sr620(1);
    end
    fopen(sr620);
    
    % Write configuration commands to SR620
    for i = 1:length(cmds)
        fprintf(sr620, cmds(i));
    end

    % Disconnect from instrument object, obj1.
    fclose(sr620);

    
    % Change default out to time + cfg_filename - .cfg + .txt later if wanted
    out_filename = 'default_out.txt';
    out_dir = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Control Software\Matlab Scripts\SR620_outs\';
    out_filename = strcat(out_dir, out_filename);
    
    % Write data out to file 
    foutid = fopen(out_filename, 'wt');
    fprintf(foutid, "");
    fclose(foutid);
    
end 
