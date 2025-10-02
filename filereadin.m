clear;
clc;
close all;

%% declare variables

%material properties
%order: Al, Br, St
density = [2810;8500;8000]; %kg/m^3
c_p = [960;380;500]; %J/kgK
k = [130;115;16.2]; %W/mK
alpha = k./(density.*c_p); %thermal diffusivity

%dimensions
x0 = 0.034925; %m
diameter = 0.0254; %m
delta_x = 0.0254; %m
cross_sec = pi * ((diameter/2)^2);

%length from x0 to the point where power supply starts
%length = [x0,(0.0762:0.0127:0.1651)]; %m

%%
folder = 'ASEN3802_HeatConduction_FA25';
a = dir(fullfile(folder));   % struct of .mat files

a(1:2) = [];

for i=1:5
    data = readmatrix(fullfile(folder, a(i).name));
    % how to get voltage and amperage from file names?
    % - options include strsplit, regex, etc.
    % ultimately, we need to use the format of each file name
    % 'material'_'volts'V_'amps'mA
    b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
    % {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'


    % v = strsplit(b{2},'V'); % volts are always in the second portion
    % ampval= strsplit(b{3},'mA'); % amps are always in the third portion
    % volts(i) = str2num(v{1}); % convert string to number (vector)
    % amps(i) = str2num(ampval{1});

    if (b{1} == "Aluminum") && (b{2} == "25V") && (b{3} == "240mA")
        length = linspace(x0,0.1651,length(data(:,2)));
        p = polyfit(length,data(:,2),3);
        y1 = polyval(p,length);
        figure;
        plot(length,y1)
        xlabel("Length along rod (m)")
        ylabel("Temp (deg C)")

    % elseif (b{1} == "Aluminum") && (b{2} == "30V") && (b{3} == "290mA")
    % elseif (b{1} == "Brass") && (b{2} == "30V") && (b{3} == "285mA")
    % elseif (b{1} == "Brass") && (b{2} == "25V") && (b{3} == "237mA")
    % elseif (b{1} == "Steel") && (b{2} == "22V") && (b{3} == "203mA")
    else
        disp("Error: Incorrect material properties")
    end

end

% figure;
% print('test','-dpng','-r300')
% saveas(a,'test2',)
% line fit slope to get H
% linearly extrapolate to get T0
