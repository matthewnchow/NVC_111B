%% Connect to instruments and data file
SR620 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 16, 'Tag', '');
if isempty(SR620)
    SR620 = gpib('NI', 0, 16);
else
    fclose(SR620);
    SR620 = SR620(1);
end
fopen(SR620);

SRS384 = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 27, 'Tag', '');
if isempty(SRS384)
    SRS384 = gpib('NI', 0, 27);
else
    fclose(SRS384);
    SRS384 = SRS384(1);
end
fopen(SRS384);

ddir = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\';
dfilename = 'test-data_5-20-19_MW10ms_1s_noratio.csv';
dfile = fopen([ddir, dfilename],'w+');

%% Communicate with instruments, sweep through desired freq range
freq_i = 2.84;
freq_f = 2.89;
N = 500;
delta_freq = (freq_f - freq_i)/N;

fprintf(dfile, '%s\r\n', 'Frequency(GHz),Ratio of Cts (Laser On/off)');
for j = 1:10
    for i = 1:N
        ghz = num2str(freq_i + delta_freq * i);
        str = ['FREQ ', ghz, 'e9'];
        fprintf(SRS384, str);
        dat = SR620_BinDump(1, SR620, 0);
        fprintf(dfile, '%s', [ghz, ',']);
        fprintf(dfile, '%f\r\n', dat);
    end
end

%% Close and clean up
fclose(dfile);
fclose(SRS384);

delete(SRS384)
delete(SR620)