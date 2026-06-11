% Function to make automatically text for inclusion in legends
function leg_text = legend_text_Ss_Sigma(Sss,Sigmas)

% Round to two significant figures 
Sss = round(Sss, 2,'significant');
Sigmas = round(Sigmas, 2,'significant');

% Make sure they are column vectors
if size(Sss,2) ~= 1;     Sss = Sss';   end
if size(Sigmas,2) ~= 1;    Sigmas = Sigmas'; end

% Convert into string array
Sss = convertCharsToStrings(cellstr(num2str(Sss)));
Sigmas = convertCharsToStrings(cellstr(num2str(Sigmas)));

% Change standard form notation to something more tex-friendly
Sss = strrep(Sss ,'e-0','\times 10^');

leg_text = strcat( '$S_s=', Sss ,'$ m$^{-1}$ ($\Sigma=', Sigmas, '$)');

end