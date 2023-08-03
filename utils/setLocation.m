function [lat_long_alt] = setLocation(location)
%setLocation sets the location of the observer
%input: location which is either a label or location vector from the list
%output: location coordinates of that label
% Tested : matlab 9.10
%     By : Saurabh Pegwal                 Jan 2022
% Example: [a,b,c]=setLocation(label)

if ischar(location)==1
    Location_label = { 'Null Island', 'REACH', 'SKA', 'HERA','PAPER' , 'GBT' , 'GMRT' , 'FAST', 'ASKAP', 'ATCA', 'MWA', 'PaST', 'LOFAR', 'CHIME','SU'};

    L_coord=[0,0,0; -30.838750, 21.374920, 1150.0000; -30.72113, 21.411128,1080; -30.721339829017378, 21.4282424879626, 1080;...
        -30.722400, 21.427800, 1080;  38.43302873568111, -79.83982982125139,807.43;...
        19.091220708385247, 74.0505174019285, 680; 25.652889, 106.856778, 5080;...
        -26.69699102820478, 116.63105557519636, 380; -30.312987696087397, 149.56440104951383, 180;...
        -26.702994018015414, 116.67038222390741, 1060.00; 42.924200, 86.716000, 2500;...
        52.90156227803543, 6.849166859086357, 140; 49.320833, -119.623611, 545; -33.928323, 18.866277, 150;];

    %assert(ismember(location,Location_label), 'Unkown location. See Location_label for allowable names');
    assert(ismember(location,Location_label),'%s\n %s\n',  'Unkown location. Allowable names' ,Location_label{:});
    for i=1:length(L_coord(:,1))
        if strcmp(Location_label{i}, location)
            lat=L_coord(i,1);long=L_coord(i,2);alt=L_coord(i,3);
            lat_long_alt = [lat, long, alt];
        end
    end
else
    lat_long_alt = [location(1), location(2), location(3)];
end
end