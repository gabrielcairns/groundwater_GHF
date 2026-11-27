% Function to make automatically text for inclusion in legends
function leg_text = legend_text_K_Pe(Ks,Pes)

% Round to two significant figures 
Ks = round(Ks, 2,'significant');
Pes = round(Pes, 2,'significant');

% Make sure they are column vectors
if size(Ks,2) ~= 1;     Ks = Ks';   end
if size(Pes,2) ~= 1;    Pes = Pes'; end

if any( isbetween(Ks,0,0.001,"open") | (Ks>1000))
    % If we're dealing with small or large numbers, use exponential form
    Kchars = num2str(Ks,'%.1e');
    % But not for zero:
    if any(Ks==0)
        Kchars(Ks==0,:) = '      0';
    end
else
    % If numbers are 'O(1)', just use default settings
    Kchars = num2str(Ks);
end

if any( isbetween(Pes,0,0.001,"open") | (Pes>1000))
    % If we're dealing with small or large numbers, use exponential form
    Pechars = num2str(Pes,'%.1e');
    % But not for zero:
    if any(Pes==0)
        Pechars(Pes==0,:) = '      0';
    end
else
    % If numbers are 'O(1)', just use default settings
    Pechars = num2str(Pes);
end

% Convert into string array
Ks = convertCharsToStrings(cellstr(num2str(Kchars)));
Pes = convertCharsToStrings(cellstr(num2str(Pechars)));

% Change standard form notation to something more tex-friendly
for i = 1:length(Ks)
    Ks{i} = float_to_tex(Ks{i});
    Pes{i} = float_to_tex(Pes{i});
end


leg_text = strcat( '$K_{sb}=', Ks ,'$ ms$^{-1}$ (Pe$=', Pes, '$)');

end