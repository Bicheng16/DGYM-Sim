function output = createReport(input)

    % Store all output file names in one cell
    output = cell(3, 1);

    % Get input file name (name.dek)
    tmp = strsplit(input, '\');
    iname = char(tmp(end));    
    % Get output file name 
    iiname = strsplit(iname, '.');
        
    output{1} = char( strcat(iiname(1), '.rpt') );  % Fluid in place
    output{2} = char( strcat(iiname(1), '_prod.rsm') );  % Time step production / pressure
    output{3} = char( strcat(iiname(1), '_inj.rsm') );  % Time step production / pressure
    
end