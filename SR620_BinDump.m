function data = SR620_BinDump(samples, handle, ratio_on) %Add instrument handle
    switch nargin
        case 1
            ratio_on = 0; 
            handle = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');
        case 0
            ratio_on = 0;
            samples = 1;
            handle = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');
    end
    
    % Create the GPIB object if it does not exist otherwise use the object that was found.
    if isempty(handle)
        handle = gpib('NI', 0, 16);
    else
        fclose(handle);
        handle = handle(1);
    end
    handle.InputBufferSize = max((128 + samples * 64), handle.InputBufferSize); %Room for 2 * number of samples * 64 (samples returned as 64 bit value)
    fopen(handle);

    mode = str2num(query(handle, 'MODE?')); %Ask sr620 what mode it's in
    x1k = logical(str2num(query(handle, 'EXPD?'))); %Ask sr620 if x1000 is on
    fwrite(handle, sprintf('BDMP%d', samples));
    rawdata = double(typecast(fread(handle, samples, 'double'), 'int64')); %convert data to useable type
    fprintf(handle, 'MODE%d', mode); %Return display to showing values
    fclose(handle);

    %Convert raw number using spec'd constant multiple
    const = 1.0;
    switch mode
        case {0, 1, 2, 4}
            const = 1.05963812934 * (10^-14);
        case 3
            const =  1.24900090270331*10^(-9);
        case 5
            const =  8.3819032 * 10^(-8);
        case 6
            const = 1/256.0;
    end

    if ratio_on 
        const = const * 1/2^40;
    end
    if x1k && (mode == 3 || mode == 4)
        const = const / 1000;
    end 

    data = rawdata * const;
end
