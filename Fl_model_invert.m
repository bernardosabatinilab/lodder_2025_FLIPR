function [xF, xT]=Fl_model_invert(params, F_data, T_data)
    R=params(1);
    f_offset=params(2); 
    t1=params(3);
    t2=params(4);
    tb=params(5);
    fb=params(6);
    s=params(7);
    
    % calculate the fraction in state 1 from the fluorescence intensity
    xF=[];
    if ~isempty(F_data)
        xF=(F_data-f_offset-fb-R*s)/(s*(1-R));
    end

    % calculate the fraction in state 1 from the fluorescence lifetime
    xT=[];
    if ~isempty(T_data)
        xT=(fb*tb/s + t2*R - T_data*(fb/s+R))./(T_data-t1+R*(t2-T_data));
    end

