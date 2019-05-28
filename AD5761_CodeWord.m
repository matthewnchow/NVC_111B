% Either returns three words or one, based on whether or not configdata is
% null. If configdata != null => returns reset word, config word, data word
% Else returns data word based on the range 

function words = AD5761R_CodeWord(vout, range, configdata) %For now range is assumed +/-10 V
    vref = 2.5;
    C = 4;
    M = 8;
    reset = ['00000111' '00000000' '00000000'];
    defaultConfig = ['00000100' '00000110' '11111000']; % +/-10 V range, 2's compliment
    switch nargin
        case 3
            % deliver 
            dataWord = ['00000011' int16(((vout / vref) +  C)/M * 2^16)];
            words = {reset, defaultConfig, dataWord}
        case 2
            
        case 1 % Assume the +/-10V range
            words = {['00000011' int16(((vout / vref) +  C)/M * 2^16)]}
    end
end 