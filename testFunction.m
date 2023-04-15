%% Initialization
close all; clearvars; clc

[filename, pathname] = uigetfile({'*.dat'},'Select APC propeller performance file');
apcPropPerfData = readAPCperf([pathname,filename]);

%% Plot data
dataset = size(apcPropPerfData,2);

rpm = zeros(dataset,1);
for idx = 1:dataset
    rpm(idx) = unique(apcPropPerfData{idx}.RPM);
end


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.J,apcPropPerfData{idx}.CT,'LineWidth',2)
end
hold off, grid on
xlabel('J'), ylabel('CT'), title('Thrust coefficient')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.J,apcPropPerfData{idx}.CP,'LineWidth',2)
end
hold off, grid on
xlabel('J'), ylabel('CP'), title('Power coefficient')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.J,apcPropPerfData{idx}.eff,'LineWidth',2)
end
hold off, grid on
xlabel('J'), ylabel('\eta'), title('Propeller efficiency')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.V_kph,apcPropPerfData{idx}.T_N,'LineWidth',2)
end
hold off, grid on
xlabel('V (km/h)'), ylabel('T (N)'), title('Propeller thrust')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.V_kph,apcPropPerfData{idx}.Q_Nm,'LineWidth',2)
end
hold off, grid on
xlabel('V (km/h)'), ylabel('Q (Nm)'), title('Propeller torque')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.V_kph,apcPropPerfData{idx}.P_W,'LineWidth',2)
end
hold off, grid on
xlabel('V (km/h)'), ylabel('P (W)'), title('Propeller power')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})

%% Data interpolation initialization

% Concatenate all datasets vertically and add RPM as final column
allDataset = [];
for idx = 1:dataset
    temp = table2array(apcPropPerfData{idx});
    allDataset = [allDataset; temp];
end

% Enable a correspondence between variable names and their positions
header = apcPropPerfData{1}.Properties.VariableNames;

% Use headFun function to link numeric array position to table variable 
% names (and hopefully improve code readability and avoid errors)

%% Interpolate among curves at given RPM
queryRPM = 6550; % desired RPM performance
queryVelocity = linspace(0,max(allDataset(:,headFun(header,'V_kph')))); % query velocity array in km/h
% queryAdvRatio = linspace(0,max(allDataset(:,headFun(header,'J')))); % query advance ratio J array

% Thrust interpolant
fThrust = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'T_N')));

% Interpolate thrust among velocity and RPM
thrust = fThrust(queryVelocity,repmat(queryRPM,1,100));
posDataIndex = thrust > 0; % do not show data with negative thrust
interpThrust = [queryVelocity(posDataIndex)', thrust(posDataIndex)'];

% Plot interpolated data
figure
plot(interpThrust(:,1), interpThrust(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Velocity (km/h)'), ylabel('Thrust (N)')
title(['Interpolated thrust at ', num2str(queryRPM), ' RPM'])


% Torque interpolant
fTorque = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'Q_Nm')));

% Interpolate torque among velocity and RPM
torque = fTorque(queryVelocity,repmat(queryRPM,1,100));
posDataIndex = torque > 0; % do not show data with negative torque
interpTorque = [queryVelocity(posDataIndex)', torque(posDataIndex)'];

% Plot interpolated data
figure
plot(interpTorque(:,1), interpTorque(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Velocity (km/h)'), ylabel('Torque (Nm)')
title(['Interpolated torque at ', num2str(queryRPM), ' RPM'])


% Power interpolant
fPower = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'P_W')));

% Interpolate power among velocity and RPM
power = fPower(queryVelocity,repmat(queryRPM,1,100));
posDataIndex = power > 0; % do not show data with negative torque
interpPower = [queryVelocity(posDataIndex)', power(posDataIndex)'];

% Plot interpolated data
figure
plot(interpPower(:,1), interpPower(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Velocity (km/h)'), ylabel('Power (W)')
title(['Interpolated power at ', num2str(queryRPM), ' RPM'])


%% Interpolate among curves at given airspeed
queryRPM = linspace(0,max(allDataset(:,headFun(header,'RPM')))); % desired RPM array
queryVelocity = 70; % query desired velocity in km/h
% queryAdvRatio = linspace(0,max(allDataset(:,headFun(header,'J')))); % query advance ratio J array

% Thrust interpolant
fThrust = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'T_N')));

% Interpolate thrust among velocity and RPM
thrust = fThrust(repmat(queryVelocity,1,100),queryRPM);
posDataIndex = thrust > 0; % do not show data with negative thrust
interpThrust = [queryRPM(posDataIndex)', thrust(posDataIndex)'];

% Plot interpolated data
figure
plot(interpThrust(:,1), interpThrust(:,2), 'k', 'LineWidth',2)
grid on, xlabel('RPM'), ylabel('Thrust (N)')
title(['Interpolated thrust at ', num2str(queryVelocity), ' km/h'])


% Torque interpolant
fTorque = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'Q_Nm')));

% Interpolate torque among velocity and RPM
torque = fTorque(repmat(queryVelocity,1,100),queryRPM);
posDataIndex = torque > 0; % do not show data with negative thrust
interpTorque = [queryRPM(posDataIndex)', torque(posDataIndex)'];

% Plot interpolated data
figure
plot(interpTorque(:,1), interpTorque(:,2), 'k', 'LineWidth',2)
grid on, xlabel('RPM'), ylabel('Torque (Nm)')
title(['Interpolated torque at ', num2str(queryVelocity), ' km/h'])


% Power interpolant
fPower = scatteredInterpolant(allDataset(:,headFun(header,'V_kph')), ...
    allDataset(:,headFun(header,'RPM')), ...
    allDataset(:,headFun(header,'P_W')));

% Interpolate power among velocity and RPM
power = fPower(repmat(queryVelocity,1,100),queryRPM);
posDataIndex = power > 0; % do not show data with negative thrust
interpPower = [queryRPM(posDataIndex)', power(posDataIndex)'];

% Plot interpolated data
figure
plot(interpPower(:,1), interpPower(:,2), 'k', 'LineWidth',2)
grid on, xlabel('RPM'), ylabel('Power (W)')
title(['Interpolated power at ', num2str(queryVelocity), ' km/h'])