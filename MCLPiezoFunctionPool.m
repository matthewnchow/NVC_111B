function MCLPiezoFunctionPool(varargin)
% Built with same structure as PIPiezoFunctionPool, written by Matthew
% Chow
    switch varargin{1}
        case 'Init'                            
            Init(); 
        case 'CheckErrors'
            CheckErrors()
        case 'SVO_Set'
            SVO_Set(varargin{2})
        case 'Move_Rel'
            Move_Rel(varargin{2})
        case 'Move_Abs'
            Move_Abs(varargin{2})
        case 'Real_Position_Q'
            Real_Pos_Q()
        case 'Curr_Voltage_Q'
            Curr_V_Q()
    end
end


function Init()
    global MCLPiezo
    
    

end
