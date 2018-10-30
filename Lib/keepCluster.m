function OUT = keepCluster
    choice = questdlg('Keep the cluster?', ...
        'Cluster Keep Menu', ...
        'Yes','No','Skip', 'Skip');
    % Handle response
    switch choice
        case 'Yes'
            disp(choice)
            OUT = 1;
        case 'No'
            disp(choice)
            OUT = 0;
        case 'Skip'
            disp('All the rest are rejected.')
            OUT = 2;
    end
end