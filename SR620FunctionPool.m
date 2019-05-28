function varargout = SR620FunctionPool(varargin)

    switch varargin{1}
        case 'GetCounts' % Works 3/30/19
            switch length(varargin)
                case 2
                    varargout = GetCounts(varargin{2});
                case 1
                    varargout = GetCounts();
            end

        case 'GetTimes' % Works 3/30/19
            switch length(varargin)
                case 2
                    varargout = GetTimes(varargin{2});
                case 1
                    varargout = GetTimes();
            end
        case 'ObjHandle' % Works 3/30/19
            varargout = {ObjHandle()};
        case 'CloseCounter' % Works 3/30/19
            fclose(ObjHandle());
        case 'ConfigFromFile' % Works 3/30/19
            switch length(varargin)
                case 3
                    ConfigFromFile(varargin{2}, varargin{3});
                case 2
                    ConfigFromFile(varargin{2});
                case 1
                    ConfigFromFile();
            end
        case 'SetGatedCounter' % Works 3/30/19
            ConfigFromFile('default_count.cfg')
    end
    
    function counts = GetCounts(size) %Sampling time   
        switch nargin
            case 0
                size = 1;
        end
        fprintf(ObjHandle(), 'MODE6');
        counts = BinDump(size);
    end 
    

    function ts = GetTimes(size)
        switch nargin
            case 0
                size = 1000;
        end
        fprintf(ObjHandle(), 'MODE0');
        ts = BinDump(size);
    end


    function data = BinDump(samples) %Make handle more global
        switch nargin
            case 0
                samples = 1;
        end
        
        sr620 = ObjHandle();

        mode = str2num(query(sr620, 'MODE?')); %#ok<*ST2NM> %Ask sr620 what mode it's in
        x1k = false;

        fwrite(sr620, sprintf('BDMP%d', samples));
        rawdata = double(typecast(fread(sr620, samples, 'double'), 'int64')); %convert data to useable type
        fprintf(sr620, 'MODE%d', mode); %Return display to showing values
        
        %Convert raw number using spec'd constant multiple
        const = 1.05963812934 * (10^-14);
        switch mode
            case 6
                const = 1/256;
            case 5
                const =  8.3819032 * 10^(-8);
            case 4
                x1k = logical(str2num(query(sr620, 'EXPD?'))); %Ask sr620 if x1000 is on
            case 3
                const =  1.24900090270331*10^(-9);
                x1k = logical(str2num(query(sr620, 'EXPD?'))); %Ask sr620 if x1000 is on
        end

        if x1k && (mode == 3 || mode == 4)
            const = const / 1000;
        end 

        data = {rawdata * const};
    end

    function cmds = readcmds(pathtofile)
    % Read configuration file and prepare list of commands
        fid = fopen(pathtofile);
        tline = string(fgetl(fid)); %Header line
        i = 1;
        cmds = [];
        indices = [1];
        header = 'Output for: ';
        firstline = true; 
        while ~strcmp(tline, "-1")
            if firstline
                header = strcat(header, tline, '\r\n');
                firstline = false;
            else
                cmds = [cmds tline];
                i = i + length(tline);
                indices = [indices, i];
            end
            tline = string(fgetl(fid));
        end
        fclose(fid);
    end 


    function ConfigFromFile(cfg_filename, cfg_dir)
        switch nargin
            case 0
                cfg_filename = 'default_count.cfg';
                cfg_dir = strcat(pwd, '\SR620_configs');
            case 1
                cfg_dir = strcat(pwd, '\SR620_configs');
        end
        
        pathtofile = strcat(cfg_dir, '\', cfg_filename);
        cmds = readcmds(pathtofile);
        
        sr620 = ObjHandle();

        % Write configuration commands to SR620
        for i = 1:length(cmds)
            fprintf(sr620, cmds(i));
        end

    end


    function hand = ObjHandle() %Want to make a global variable and update it eventually
        % Find a GPIB object. Call it sr620
        sr620 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');

        % Create the GPIB object if it does not exist
        % otherwise use the object that was found.
        if isempty(sr620)
            sr620 = gpib('NI', 0, 16);
        else
            fclose(sr620);
            sr620 = sr620(1);
        end
        sr620.InputBufferSize = 8 * 4 * 10^6; %Room for 4 Mega measurements
        fopen(sr620);
        hand = sr620;
    end

end 
