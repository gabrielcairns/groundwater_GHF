function str = float_to_tex(str)
% Function to convert string numbers in standard form from floating point
% into TeX style
% Fails if string exponents are 3 digits

    exppos = strfind(str,'e-');
    if ~isempty(exppos)
        if str(exppos+2) == '0'
            exponent = str(exppos+3);
            str = insertAfter(str,exppos+3,'}');
            str = strrep(str ,'e-0','\times 10^{-');
        else
            exponent = str([exppos+2 exppos+3]);
            str = insertAfter(str,exppos+3,'}');
            str = strrep(str ,'e-','\times 10^{-');
        end
        
    end

    exppos = strfind(str,'e+');
    if ~isempty(exppos)
        if str(exppos+2) == 0
            exponent = str(exppos+3);
            str = insertAfter(str,exppos+3,'}');
            str = strrep(str ,'e+0','\times 10^{-');
        else
            exponent = str([exppos+2 exppos+3]);
            str = insertAfter(str,exppos+3,'}');
            str = strrep(str ,'e+','\times 10^{-');
        end

    end