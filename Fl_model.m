function [F, T]=Fl_model(params, f)
    R=params(1);
    f_offset=params(2); 
    t1=params(3);
    t2=params(4);
    tb=params(5);
    fb=params(6);
    s=params(7);

    % calculate the fluorescence intensity
    F = f_offset + fb + s * (f + R*(1-f)); 

    % calculate the fluorescence lifetime
    T = (fb*tb/s + t1*f + t2*R*(1-f)) ./ (fb/s + f + R*(1-f)); 
