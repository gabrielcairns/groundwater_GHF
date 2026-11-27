% Function to make automatically text for inclusion in legends
function leg_text = legend_text_Ss_Sigma(Sss,Sigmas)

% Round to two significant figures 
Sss = round(Sss, 2,'significant');
Sigmas = round(Sigmas, 2,'significant');

% Make sure they are column vectors
if size(Sss,2) ~= 1;     Sss = Sss';   end
if size(Sigmas,2) ~= 1;    Sigmas = Sigmas'; end

if any( isbetween(Sss,0,0.001,"open") | (Sss>1000))
    % If we're dealing with small or large numbers, use exponential form
    Sschars = num2str(Sss,'%.1e');
    % But not for zero:
    if any(Sss==0)
        Sschars(Sss==0,:) = '      0';
    end
else
    % If numbers are 'O(1)', just use default settings
    Sschars = num2str(Sss);
end

if any( isbetween(Sigmas,0,0.001,"open") | (Sigmas>1000))
    % If we're dealing with small or large numbers, use exponential form
    Sigmachars = num2str(Sigmas,'%.1e');
    % But not for zero:
    if any(Sigmas==0)
        Sigmachars(Sigmas==0,:) = '      0';
    end
else
    % If numbers are 'O(1)', just use default settings
    Sigmachars = num2str(Sigmas);
end


% Convert into string array
Sss = convertCharsToStrings(cellstr(num2str(Sschars)));
Sigmas = convertCharsToStrings(cellstr(num2str(Sigmachars)));

% Change standard form notation to something more tex-friendly
for i = 1:length(Sss)
    Sss{i} = float_to_tex(Sss{i});
    Sigmas{i} = float_to_tex(Sigmas{i});
end
leg_text = strcat( '$S_s=', Sss ,'$ m$^{-1}$ ($\Sigma=', Sigmas, '$)');

end