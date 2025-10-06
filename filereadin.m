
clear;
clc;
close all;

%Author: Parinie Gupta

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
delta_x = 0.0127; %m
cross_sec = pi * ((diameter/2)^2);
rod_length = 0.180975; % m
%thermocouple locations
tc_loc = 0.0762 + (0:7)*delta_x;

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
        %Task 1
        %NOTE: Polyfit gives 2 arguments: first is slope and second is
        %intercept

        %defining entire length of rod
        l_steady = linspace(0,rod_length,100);

        %steady state temp
        steady_temp = [data(end-2,2:9)]+273.15; %K
        
        %fitting
        [p1,S1] = polyfit(tc_loc,steady_temp,1);
        [y1,delta1] = polyval(p1,l_steady,S1);

        %calculating H_exp
        H_exp(i) = p1(1); %K/m

        %extract T0 from polyfit (intercept)
        T0(i) = p1(2); %K

        %calculating H_an
        %H_an = VI/kA
        H_an(i) = ((volts(i) * (amps(i)/1000)) / (k(1) * cross_sec)); %K/m

        %plotting steady state temp starting from x0
        figure;
        hold on;
        grid("on")
        %plotting with H_exp slope
        plot(l_steady,y1,LineWidth=2)
        plot(l_steady,y1+delta1,'r--',LineWidth=2)
        plot(l_steady,y1-delta1,'r--',LineWidth=2)
        %plot with H_an slope
        plot(l_steady,(H_an(i).*l_steady)+T0(i),'k',LineWidth=2)
        %marking TC locations
        plot(tc_loc, steady_temp, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('H_{exp}','+error','-error','H_{an}','TC Data')
        title('Aluminum 25V 240mA')
        subtitle('Steady State')

        %saving figures
        fname = sprintf('steady_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

        %Task 2
        %assuming T0 is the first TC temp
        T0_init = data(1,2)+273.15;%K

        %temp from thermocouples at t=0
        temp_init = [data(1,2:9)]+273.15;%K

        %defining entire length of rod
        l_init = linspace(0,rod_length,100);

        %fitting
        [p2,S2] = polyfit(tc_loc,temp_init,1);
        [y2,delta2] = polyval(p2,l_init,S2);

        %extracting M_exp
        M_exp(i) = p2(1);%K/m

        %plot
        figure;
        grid("on")
        hold on;
        %plotting with M_exp slope
        plot(l_init,y2,LineWidth=2)
        plot(l_init,y2+delta2,'r--',LineWidth=2)
        plot(l_init,y2-delta2,'r--',LineWidth=2)
        %marking TC locations
        plot(tc_loc, temp_init, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        ylim([288,292])
        legend('M_{exp}','+error','-error','TC Data')
        title('Aluminum 25V 240mA')
        subtitle('Initial State')

        %saving figures
        fname = sprintf('initial_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

    elseif (b{1} == "Aluminum") && (b{2} == "30V") && (b{3} == "290mA")
        %Task 1
        %NOTE: Polyfit gives 2 arguments: first is slope and second is
        %intercept

        %defining entire length of rod
        l_steady = linspace(0,rod_length,100);

        %steady state temp
        steady_temp = [data(end-2,2:9)]+273.15; %K
        
        %fitting
        [p1,S1] = polyfit(tc_loc,steady_temp,1);
        [y1,delta1] = polyval(p1,l_steady,S1);

        %calculating H_exp
        H_exp(i) = p1(1); %K/m

        %extract T0 from polyfit (intercept)
        T0(i) = p1(2); %K

        %calculating H_an
        %H_an = VI/kA
        H_an(i) = ((volts(i) * (amps(i)/1000)) / (k(1) * cross_sec)); %K/m

        %plotting steady state temp starting from x0
        figure;
        hold on;
        grid("on")
        %plotting with H_exp slope
        plot(l_steady,y1,LineWidth=2)
        plot(l_steady,y1+delta1,'r--',LineWidth=2)
        plot(l_steady,y1-delta1,'r--',LineWidth=2)
        %plot with H_an slope
        plot(l_steady,(H_an(i).*l_steady)+T0(i),'k',LineWidth=2)
        %marking TC locations
        plot(tc_loc, steady_temp, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('H_{exp}','+error','-error','H_{an}','TC Data')
        title('Aluminum 30V 290mA')
        subtitle('Steady State')

        %saving figures
        fname = sprintf('steady_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

        %Task 2
        %assuming T0 is the first TC temp
        T0_init = data(1,2)+273.15;%K

        %temp from thermocouples at t=0
        temp_init = [data(1,2:9)]+273.15;%K

        %defining entire length of rod
        l_init = linspace(0,rod_length,100);

        %fitting
        [p2,S2] = polyfit(tc_loc,temp_init,1);
        [y2,delta2] = polyval(p2,l_init,S2);

        %extracting M_exp
        M_exp(i) = p2(1);%K/m

        %plot
        figure;
        grid("on")
        hold on;
        %plotting with M_exp slope
        plot(l_init,y2,LineWidth=2)
        plot(l_init,y2+delta2,'r--',LineWidth=2)
        plot(l_init,y2-delta2,'r--',LineWidth=2)
        %marking TC locations
        plot(tc_loc, temp_init, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        ylim([288,292])
        legend('M_{exp}','+error','-error','TC Data')
        title('Aluminum 30V 290mA')
        subtitle('Initial State')

        %saving figures
        fname = sprintf('initial_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

    elseif (b{1} == "Brass") && (b{2} == "25V") && (b{3} == "237mA")
        %Task 1
        %NOTE: Polyfit gives 2 arguments: first is slope and second is
        %intercept

        %defining entire length of rod
        l_steady = linspace(0,rod_length,100);

        %steady state temp
        steady_temp = [data(end-2,2:9)]+273.15; %K
        
        %fitting
        [p1,S1] = polyfit(tc_loc,steady_temp,1);
        [y1,delta1] = polyval(p1,l_steady,S1);

        %calculating H_exp
        H_exp(i) = p1(1); %K/m

        %extract T0 from polyfit (intercept)
        T0(i) = p1(2); %K

        %calculating H_an
        %H_an = VI/kA
        H_an(i) = ((volts(i) * (amps(i)/1000)) / (k(2) * cross_sec)); %K/m

        %plotting steady state temp starting from x0
        figure;
        hold on;
        grid("on")
        %plotting with H_exp slope
        plot(l_steady,y1,LineWidth=2)
        plot(l_steady,y1+delta1,'r--',LineWidth=2)
        plot(l_steady,y1-delta1,'r--',LineWidth=2)
        %plot with H_an slope
        plot(l_steady,(H_an(i).*l_steady)+T0(i),'k',LineWidth=2)
        %marking TC locations
        plot(tc_loc, steady_temp, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('H_{exp}','+error','-error','H_{an}','TC Data')
        title('Brass 25V 237mA')
        subtitle('Steady State')

        %saving figures
        fname = sprintf('steady_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

        %Task 2
        %assuming T0 is the first TC temp
        T0_init = data(1,2)+273.15;%K

        %temp from thermocouples at t=0
        temp_init = [data(1,2:9)]+273.15;%K

        %defining entire length of rod
        l_init = linspace(0,rod_length,100);

        %fitting
        [p2,S2] = polyfit(tc_loc,temp_init,1);
        [y2,delta2] = polyval(p2,l_init,S2);

        %extracting M_exp
        M_exp(i) = p2(1);%K/m

        %plot
        figure;
        grid("on")
        hold on;
        %plotting with M_exp slope
        plot(l_init,y2,LineWidth=2)
        plot(l_init,y2+delta2,'r--',LineWidth=2)
        plot(l_init,y2-delta2,'r--',LineWidth=2)
        %marking TC locations
        plot(tc_loc, temp_init, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('M_{exp}','+error','-error','TC Data')
        title('Brass 25V 237mA')
        subtitle('Initial State')

        %saving figures
        fname = sprintf('initial_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')
        
    elseif (b{1} == "Brass") && (b{2} == "30V") && (b{3} == "285mA")
        %Task 1
        %NOTE: Polyfit gives 2 arguments: first is slope and second is
        %intercept

        %defining entire length of rod
        l_steady = linspace(0,rod_length,100);

        %steady state temp
        steady_temp = [data(end-2,2:9)]+273.15; %K
        
        %fitting
        [p1,S1] = polyfit(tc_loc,steady_temp,1);
        [y1,delta1] = polyval(p1,l_steady,S1);

        %calculating H_exp
        H_exp(i) = p1(1); %K/m

        %extract T0 from polyfit (intercept)
        T0(i) = p1(2); %K

        %calculating H_an
        %H_an = VI/kA
        H_an(i) = ((volts(i) * (amps(i)/1000)) / (k(2) * cross_sec)); %K/m

        %plotting steady state temp starting from x0
        figure;
        hold on;
        grid("on")
        %plotting with H_exp slope
        plot(l_steady,y1,LineWidth=2)
        plot(l_steady,y1+delta1,'r--',LineWidth=2)
        plot(l_steady,y1-delta1,'r--',LineWidth=2)
        %plot with H_an slope
        plot(l_steady,(H_an(i).*l_steady)+T0(i),'k',LineWidth=2)
        %marking TC locations
        plot(tc_loc, steady_temp, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('H_{exp}','+error','-error','H_{an}','TC Data')
        title('Brass 30V 285mA')
        subtitle('Steady State')

        %saving figures
        fname = sprintf('steady_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

        %Task 2
        %assuming T0 is the first TC temp
        T0_init = data(1,2)+273.15;%K

        %temp from thermocouples at t=0
        temp_init = [data(1,2:9)]+273.15;%K

        %defining entire length of rod
        l_init = linspace(0,rod_length,100);

        %fitting
        [p2,S2] = polyfit(tc_loc,temp_init,1);
        [y2,delta2] = polyval(p2,l_init,S2);

        %extracting M_exp
        M_exp(i) = p2(1);%K/m

        %plot
        figure;
        grid("on")
        hold on;
        %plotting with M_exp slope
        plot(l_init,y2,LineWidth=2)
        plot(l_init,y2+delta2,'r--',LineWidth=2)
        plot(l_init,y2-delta2,'r--',LineWidth=2)
        %marking TC locations
        plot(tc_loc, temp_init, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('M_{exp}','+error','-error','TC Data')
        title('Brass 30V 285mA')
        subtitle('Initial State')

        %saving figures
        fname = sprintf('initial_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

    elseif (b{1} == "Steel") && (b{2} == "22V") && (b{3} == "203mA")
        %Task 1
        %NOTE: Polyfit gives 2 arguments: first is slope and second is
        %intercept

        %defining entire length of rod
        l_steady = linspace(0,rod_length,100);

        %steady state temp
        steady_temp = [data(end-2,2:9)]+273.15; %K
        
        %fitting
        [p1,S1] = polyfit(tc_loc,steady_temp,1);
        [y1,delta1] = polyval(p1,l_steady,S1);

        %calculating H_exp
        H_exp(i) = p1(1); %K/m

        %extract T0 from polyfit (intercept)
        T0(i) = p1(2); %K

        %calculating H_an
        %H_an = VI/kA
        H_an(i) = ((volts(i) * (amps(i)/1000)) / (k(3) * cross_sec)); %K/m

        %plotting steady state temp starting from x0
        figure;
        hold on;
        grid("on")
        %plotting with H_exp slope
        plot(l_steady,y1,LineWidth=2)
        plot(l_steady,y1+delta1,'r--',LineWidth=2)
        plot(l_steady,y1-delta1,'r--',LineWidth=2)
        %plot with H_an slope
        plot(l_steady,(H_an(i).*l_steady)+T0(i),'k',LineWidth=2)
        %marking TC locations
        plot(tc_loc, steady_temp, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('H_{exp}','+error','-error','H_{an}','TC Data')
        title('Steel 22V 203mA')
        subtitle('Steady State')

        %saving figures
        fname = sprintf('steady_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')

        %Task 2
        %assuming T0 is the first TC temp
        T0_init = data(1,2)+273.15;%K

        %temp from thermocouples at t=0
        temp_init = [data(1,2:9)]+273.15;%K

        %defining entire length of rod
        l_init = linspace(0,rod_length,100);

        %fitting
        [p2,S2] = polyfit(tc_loc,temp_init,1);
        [y2,delta2] = polyval(p2,l_init,S2);

        %extracting M_exp
        M_exp(i) = p2(1);%K/m

        %plot
        figure;
        grid("on")
        hold on;
        %plotting with M_exp slope
        plot(l_init,y2,LineWidth=2)
        plot(l_init,y2+delta2,'r--',LineWidth=2)
        plot(l_init,y2-delta2,'r--',LineWidth=2)
        %marking TC locations
        plot(tc_loc, temp_init, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
        xlabel("Length along rod (m)")
        ylabel("Temp (K)")
        legend('M_{exp}','+error','-error','TC Data')
        title('Steel 22V 203mA')
        subtitle('Initial State')

        %saving figures
        fname = sprintf('initial_%s_%s_%s.png', b{1}, b{2}, b{3});
        saveas(gcf,fname,'png')
        
    else
        disp("Error: Incorrect material properties")
    end

end
