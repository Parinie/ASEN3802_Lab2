clear;
clc;
close all;

%% declare variables

%material properties
%order: Al, Br, St
density = [2810;8500;8000]; %kg/m^3
c_p = [960;380;500]; %J/kgK
k = [130;115;16.2]; %W/mK themal conductivity
alpha = k./(density.*c_p); %thermal diffusivity

%dimensions
x0 = 0.034925; %m
diameter = 0.0254; %m
delta_x = 0.0254; %m
cross_sec = pi * ((diameter/2)^2);


%%
folder = '3802Lab2_data';
a = dir(fullfile(folder));   % struct of .mat files

a(1:2) = [];

for i=1:5
    data = readmatrix(fullfile(folder, a(i).name));
    
    b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3

    % {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
    v = strsplit(b{2},'V'); % volts are always in the second portion
    ampval= strsplit(b{3},'mA'); % amps are always in the third portion
    volts(i) = str2num(v{1}); % convert string to number (vector)
    amps(i) = str2num(ampval{1});

    if (b{1} == "Aluminum") && (b{2} == "25V") && (b{3} == "240mA")
        %temp from thermocouples at t=0
        temp_time0 = [data(1,2),data(1,2),data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8),data(1,9)];
        length = linspace(0,0.2159,length(temp_time0));

        %fitting
        p = polyfit(length,temp_time0,1);
        y1 = polyval(p,length);

        %plot
        figure;
        hold on;
        plot(length,y1)
        xlabel("Length along rod (m)")
        ylabel("Temp (deg C)")
        title('Aluminum 25V 240mA')

        %extract T0
        T0 = y1(1); %deg C

        %calculating H_an
        %H_an = VI/kA
        H_an = (volts(i) * amps(i)) / (k(1) * cross_sec);

        %calculating H_exp
        %steady state temp
        % steady_temp = [data(end,2),data(end,2),data(end,3),data(end,4),data(end,5),data(end,6),data(end,7),data(end,8)];
        % p = polyfit(length,steady_temp,1);
        % y1 = polyval(p,length);
        % figure;
        % plot(length,y1)

        

    % elseif (b{1} == "Aluminum") && (b{2} == "30V") && (b{3} == "290mA")
    %     %temp from thermocouples at t=0
    %     temp_time0 = [data(1,2),data(1,2),data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8),data(1,9)];
    %     length = linspace(0,0.2159,length(temp_time0));
    % 
    %     %fitting
    %     p = polyfit(length,temp_time0,1);
    %     y1 = polyval(p,length);
    % 
    %     %extract T0
    %     T0 = y1(1); %deg C
    % 
    %     %plot
    %     figure;
    %     hold on;
    %     plot(length,y1)
    %     xlabel("Length along rod (m)")
    %     ylabel("Temp (deg C)")
    %     title('Aluminum 30V 290mA')
    % 
    %     %calculating H_an
    % elseif (b{1} == "Brass") && (b{2} == "30V") && (b{3} == "285mA")
    %     %temp from thermocouples at t=0
    %     temp_time0 = [data(1,2),data(1,2),data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8),data(1,9)];
    %     length = linspace(0,0.2159,length(temp_time0));
    % 
    % 
    %     %fitting
    %     p = polyfit(length,temp_time0,1);
    %     y1 = polyval(p,length);
    % 
    %     %extract T0
    %     T0 = y1(1); %deg C
    % 
    %     %plot
    %     figure;
    %     hold on;
    %     plot(length,y1)
    %     xlabel("Length along rod (m)")
    %     ylabel("Temp (deg C)")
    %     title('Brass 25V 237mA')

        %calculating H_an
    % elseif (b{1} == "Brass") && (b{2} == "25V") && (b{3} == "237mA")
    %     %temp from thermocouples at t=0
    %     temp_time0 = [data(1,2),data(1,2),data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8),data(1,9)];
    % 
    %     %fitting
    %     p = polyfit(length,temp_time0,2);
    %     y1 = polyval(p,length);
    % 
    %     %extract T0
    %     T0 = y1(1); %deg C
    % 
    %     %plot
    %     figure;
    %     hold on;
    %     plot(length,y1)
    %     xlabel("Length along rod (m)")
    %     ylabel("Temp (deg C)")
    % 
    %     %calculating H_an
    % elseif (b{1} == "Steel") && (b{2} == "22V") && (b{3} == "203mA")
    %     %temp from thermocouples at t=0
    %     temp_time0 = [data(1,2),data(1,2),data(1,3),data(1,4),data(1,5),data(1,6),data(1,7),data(1,8),data(1,9)];
    % 
    %     %fitting
    %     p = polyfit(length,temp_time0,2);
    %     y1 = polyval(p,length);
    % 
    %     %extract T0
    %     T0 = y1(1); %deg C
    % 
    %     %plot
    %     figure;
    %     hold on;
    %     plot(length,y1)
    %     xlabel("Length along rod (m)")
    %     ylabel("Temp (deg C)")
    % 
    %     %calculating H_an
    % else
    %     disp("Error: Incorrect material properties")
    end

end

% figure;
% print('test','-dpng','-r300')
% saveas(a,'test2',)
% line fit slope to get H
% linearly extrapolate to get T0
