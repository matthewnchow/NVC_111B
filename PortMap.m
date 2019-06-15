function path=PortMap(what)

switch what
    case 'BackupFile'
        path='C:\support\NVC\MATLAB_Code_111Lab\Data\TempDataBackup\Temp.mat';
    case 'BackupImageFile'
        path='C:\support\NVC\MATLAB_Code_111Lab\Data\TempDataBackup\TempImage.mat';
    case 'Ctr src'
       path='/Dev1/PFI8';
    case {'Ctr StartStop','StartStop'}
        path = '/Dev1/PFI9';
    case {'Ctr ArmStart', 'ArmStart'}
        path = '/Dev1/PFI7';
    case 'myDAQ trig'
        path = '/myDAQ1/DIO0';
    case 'Ctr sampclk'
        path  = PortMap('Ctr Trig'); %'/Dev1/PFI1'; 
    case 'Ctr gate'
        path='/Dev1/PFI9';
    case 'Ctr Trig'
        path='/Dev1/PFI13'; %Reserve PFI13 for pulse train generation ctr1? 
    case 'Ctr in'
        path='Dev1/ctr0';
    case 'Ctr out'
        path='Dev1/ctr1';
    case 'Galvo x'
        path='Dev1/ao0';
    case 'Galvo y'
        path='Dev1/ao1';
    case {'MCL Piezo', 'MCL piezo', 'Obj_Piezo'}'
        path = 'myDAQ1/ao0';
%     case 'APD in' % We don't have an APD
%         path='Dev1/ai0';
    case 'SG ext mod'
        path='myDAQ1/ao1';
    case 'SG mod output'
        path = 'Dev1/ai0';
    case 'SG com'
        path='com13';
    case 'meas'
        path='SPCM';
    case 'gConfocal path'
        path='C:\support\NVC\MATLAB_Code_111Lab\Sets\';
    case 'Data'
        path='C:\support\NVC\MATLAB_Code_111Lab\Data\';
    case 'LF save'
        path='C:\support\NVC\MATLAB_Code_111Lab\Data\LightField\';
    case {'C10', 'Cyc10', 'PB'}
        path = 'COM12'; %This is RS232 to USB com connection to FPGA
end