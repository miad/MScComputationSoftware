letters=['A' 'B' 'C' 'D' 'E' 'F'];
parfor lett=letters
    MEigen(['matrix' lett '.dat'],['eigenv' lett '.dat']);
end