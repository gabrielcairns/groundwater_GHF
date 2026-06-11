% Function to make automatically text for inclusion in legends
function leg_text = legend_text_K_Pe(Ks,Pes)

% Round to two significant figures 
Ks = round(Ks, 2,'significant');
Pes = round(Pes, 2,'significant');

% Make sure they are column vectors
if size(Ks,2) ~= 1;     Ks = Ks';   end
if size(Pes,2) ~= 1;    Pes = Pes'; end

% Convert into string array
Ks = convertCharsToStrings(cellstr(num2str(Ks)));
Pes = convertCharsToStrings(cellstr(num2str(Pes)));

% Change standard form notation to something more tex-friendly
Ks = strrep(Ks ,'e-0','\times 10^');

leg_text = strcat( '$K_{sb}=', Ks ,'$ ms$^{-1}$ (Pe$=', Pes, '$)');

end