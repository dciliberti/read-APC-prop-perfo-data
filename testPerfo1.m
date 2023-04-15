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

xVector = allDataset(:,headFun(header,'V_kph'));
yVector = allDataset(:,headFun(header,'RPM'));

[xNorm, Cx, Sx] = normalize(xVector);
[yNorm, Cy, Sy] = normalize(yVector);


% Thrust interpolant
vVector = allDataset(:,headFun(header,'T_N'));
[vNorm, Cv, Sv] = normalize(vVector);
fThrust = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate thrust among velocity and RPM
thrust = fThrust((queryVelocity-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = thrust*Sv+Cv > 0; % do not show data with negative thrust
interpThrust = [queryVelocity(posDataIndex)', thrust(posDataIndex)'*Sv+Cv];


% Torque interpolant
vVector = allDataset(:,headFun(header,'Q_Nm'));
[vNorm, Cv, Sv] = normalize(vVector);
fTorque = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate torque among velocity and RPM
torque = fTorque((queryVelocity-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = torque*Sv+Cv > 0; % do not show data with negative torque
interpTorque = [queryVelocity(posDataIndex)', torque(posDataIndex)'*Sv+Cv];


% Power interpolant
vVector = allDataset(:,headFun(header,'P_W'));
[vNorm, Cv, Sv] = normalize(vVector);
fPower = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate power among velocity and RPM
power = fPower((queryVelocity-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = power*Sv+Cv > 0; % do not show data with negative torque
interpPower = [queryVelocity(posDataIndex)', power(posDataIndex)'*Sv+Cv];


% Plot interpolated data
figure
subplot(3,1,1)
plot(interpThrust(:,1), interpThrust(:,2), 'k', 'LineWidth',2)
grid on, ylabel('Thrust (N)')

subplot(3,1,2)
plot(interpTorque(:,1), interpTorque(:,2), 'k', 'LineWidth',2)
grid on, ylabel('Torque (Nm)')

subplot(3,1,3)
plot(interpPower(:,1), interpPower(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Velocity (Km/h)'), ylabel('Power (W)')

sgtitle(['Interpolated curves at ', num2str(queryRPM), ' RPM'])


%% Interpolate among curves at given RPM
queryAdvRatio = linspace(0,max(allDataset(:,headFun(header,'J')))); % query advance ratio J array

xVector = allDataset(:,headFun(header,'J'));
yVector = allDataset(:,headFun(header,'RPM'));

[xNorm, Cx, Sx] = normalize(xVector);
[yNorm, Cy, Sy] = normalize(yVector);


% Thrust coefficient interpolant
vVector = allDataset(:,headFun(header,'CT'));
[vNorm, Cv, Sv] = normalize(vVector);
fThrustCoeff = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate thrust among advance ratio and RPM
thrustCoeff = fThrustCoeff((queryAdvRatio-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = thrustCoeff*Sv+Cv > 0; % do not show data with negative CT values
interpThrustCoeff = [queryAdvRatio(posDataIndex)', thrustCoeff(posDataIndex)'*Sv+Cv];


% Power coefficient interpolant
vVector = allDataset(:,headFun(header,'CP'));
[vNorm, Cv, Sv] = normalize(vVector);
fPowerCoeff = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate torque among velocity and RPM
powerCoeff = fPowerCoeff((queryAdvRatio-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = powerCoeff*Sv+Cv > 0; % do not show data with negative CP values
interpPowerCoeff = [queryAdvRatio(posDataIndex)', powerCoeff(posDataIndex)'*Sv+Cv];


% Efficiency interpolant
vVector = allDataset(:,headFun(header,'eff'));
[vNorm, Cv, Sv] = normalize(vVector);
fEfficiency = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate efficiency among velocity and RPM
efficiency = fEfficiency((queryAdvRatio-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = efficiency*Sv+Cv > 0; % do not show data with negative efficiency
interpEfficiency = [queryAdvRatio(posDataIndex)', efficiency(posDataIndex)'*Sv+Cv];


% Plot interpolated data
figure
subplot(3,1,1)
plot(interpThrustCoeff(:,1), interpThrustCoeff(:,2), 'k', 'LineWidth',2)
grid on, ylabel('Thrust coeff. C_T')

subplot(3,1,2)
plot(interpPowerCoeff(:,1), interpPowerCoeff(:,2), 'k', 'LineWidth',2)
grid on, ylabel('Power coeff. C_P')

subplot(3,1,3)
plot(interpEfficiency(:,1), interpEfficiency(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Advance ratio J'), ylabel('Prop. efficiency \eta')

sgtitle(['Interpolated curves at ', num2str(queryRPM), ' RPM'])

%% Interpolate propeller efficiency with airspeed

figure, hold on
for idx = 1:dataset
    plot(apcPropPerfData{idx}.V_kph,apcPropPerfData{idx}.eff,'LineWidth',2)
end
hold off, grid on
xlabel('V (km/h)'), ylabel('\eta'), title('Propeller efficiency')
legend({[repmat('RPM ',dataset,1),num2str(rpm)]})


xVector = allDataset(:,headFun(header,'V_kph'));
yVector = allDataset(:,headFun(header,'RPM'));

[xNorm, Cx, Sx] = normalize(xVector);
[yNorm, Cy, Sy] = normalize(yVector);

% Efficiency interpolant
vVector = allDataset(:,headFun(header,'eff'));
[vNorm, Cv, Sv] = normalize(vVector);
fEffVsAirspeed = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

% Interpolate efficiency among velocity and RPM
effVsAirspeed = fEffVsAirspeed((queryVelocity-Cx)/Sx,repmat((queryRPM-Cy)/Sy,1,100));
posDataIndex = effVsAirspeed*Sv+Cv > 0; % do not show data with negative efficiency
interpEffVsAirspeed = [queryVelocity(posDataIndex)', effVsAirspeed(posDataIndex)'*Sv+Cv];

figure
plot(interpEffVsAirspeed(:,1), interpEffVsAirspeed(:,2), 'k', 'LineWidth',2)
grid on, xlabel('Velocity (Km/h)'), ylabel('Prop. efficiency \eta')
title(['Interpolated curves at ', num2str(queryRPM), ' RPM'])