%% Initialization
close all; clearvars; clc

[filename, pathname] = uigetfile({'*.dat'},'Select APC propeller performance file');
apcPropPerfData = readAPCperf([pathname,filename]);

%% Data interpolation initialization
dataset = size(apcPropPerfData,2);

rpm = zeros(dataset,1);
for idx = 1:dataset
    rpm(idx) = unique(apcPropPerfData{idx}.RPM);
end

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

%% Interpolate among curves at several airspeeds
queryRPM = linspace(0,max(allDataset(:,headFun(header,'RPM')))); % desired RPM array

xVector = allDataset(:,headFun(header,'V_kph'));
yVector = allDataset(:,headFun(header,'RPM'));

[xNorm, Cx, Sx] = normalize(xVector);
[yNorm, Cy, Sy] = normalize(yVector);

allInterpThrust = [];
allInterpTorque = [];
allInterpPower = [];
velArray = 10:10:100; % airspeed array in km/h
for queryVelocity = velArray % cycle over velocity array

    % Thrust interpolant
    vVector = allDataset(:,headFun(header,'T_N'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fThrust = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate thrust among velocity and RPM
    thrust = fThrust(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = thrust*Sv+Cv > 0; % do not show data with negative thrust
    interpThrust = [queryRPM(posDataIndex)', thrust(posDataIndex)'*Sv+Cv'];


    % Torque interpolant
    vVector = allDataset(:,headFun(header,'Q_Nm'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fTorque = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate torque among velocity and RPM
    torque = fTorque(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = torque*Sv+Cv > 0; % do not show data with negative thrust
    interpTorque = [queryRPM(posDataIndex)', torque(posDataIndex)'*Sv+Cv];


    % Power interpolant
    vVector = allDataset(:,headFun(header,'P_W'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fPower = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate power among velocity and RPM
    power = fPower(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = power*Sv+Cv > 0; % do not show data with negative thrust
    interpPower = [queryRPM(posDataIndex)', power(posDataIndex)'*Sv+Cv];


    % Plot interpolated data
    figure(1)
    subplot(3,1,1), hold on
    plot(interpThrust(:,1), interpThrust(:,2), 'LineWidth',2)

    subplot(3,1,2), hold on
    plot(interpTorque(:,1), interpTorque(:,2), 'LineWidth',2)

    subplot(3,1,3), hold on
    plot(interpPower(:,1), interpPower(:,2), 'LineWidth',2)


    % Keep all interpolated data
    allInterpThrust = [allInterpThrust; 
                       repmat(queryVelocity,length(interpThrust),1), interpThrust];
    allInterpTorque = [allInterpTorque; 
                       repmat(queryVelocity,length(interpTorque),1), interpTorque];
    allInterpPower = [allInterpPower; 
                       repmat(queryVelocity,length(interpPower),1), interpPower];
end
figure(1)
subplot(3,1,1), hold off, grid on, ylabel('Thrust (N)')
subplot(3,1,2), hold off, grid on, ylabel('Torque (Nm)')
subplot(3,1,3), hold off, grid on, xlabel('RPM'), ylabel('Power (W)')
legend({[repmat('V (km/h) = ',length(velArray),1),num2str(velArray')]})