clear;clc;warning('off','all');rng('default');

%% Path configuration %%

path = pwd;
data_path = [erase(path,'\code') '\data'];

load([data_path '\best_individual_40_withReservoir_VIRTUAL20SENSORS.mat']);
vsensors = sensors;

load([data_path '\leaktionary_Modena_45to65_1p100_testing.mat']);

load([data_path '\ModenaNewInfo.mat']); 
load([data_path '\ModenaNewPipeDistance.mat']);
load([data_path '\ModenaNewIncidenceMatrix.mat']);
load([data_path '\best_individual_20_withReservoir.mat']);
gI = graphInfo;
G = gI.G; nc = gI.nc; A = gI.A; WA = gI.WA;
reservoirsID = gI.reservoirsID; N=gI.N; d = gI.d;
D = diag(sum(A));

precision = 100; % it translates into a precision of 0.01 
scenario_type = 'leaky';%'nominal';%

%% Running parameters %%

leak = 88; % Leak scenario (node where the leak occurs)
nom = 272; 

time = 10; % Time instant

idk = vsensors; % Demand sensor indices

t = 10.674*d.LinkLength./((d.LinkDiameter/1000).^4.871.*d.LinkRoughnessCoeff.^1.852);
T = diag(t); invT = inv(T);

% GSI-related parameters %

tau=1000;
PA = PipeDistance.*A;
W = 1./PA;
W(W==Inf) = 0;

WDeg = diag(sum(W)); WDeg1 = diag(sum(W).^-1);
WLap = WDeg - W;
WDeg2 = WDeg.^-2; WDeg2(isinf(WDeg2))=0;
WLap3 = WLap*WDeg2*WLap;

s = zeros(N,1); s(sensors) = 1; 
S = eye(N).*(s*s');

demand = Leaktionary{leak}.Demand(time,:)'; 
demand = fix(demand*precision)/precision;
demand = demand/1000; % SI - m3/s

h = Leaktionary{leak}.Head(time,:)';
q = Leaktionary{leak}.Flows(time,:)'/1000;
E = d.NodesConnectingLinksIndex; 

demand_nom = Leaktionary{nom}.Demand(time,:)'; 
demand_nom = fix(demand_nom*precision)/precision;
demand_nom = demand_nom/1000; % SI - m3/s
h_nom = Leaktionary{nom}.Head(time,:)';
q_nom = Leaktionary{nom}.Flows(time,:)'/1000;

alfa = length(idk)/length(demand);

% Precision %

demand_s = demand(idk);
h_s = fix(h(sensors)*precision)/precision;

demand_nom_s = demand_nom(idk);
h_nom_s = fix(h_nom(sensors)*precision)/precision;

%% GSI %%

hGSI = GSI(h_s,WLap3,S,I,...
    reservoirsID,...
    sensors,tau);

hGSI_nom = GSI(h_nom_s,WLap3,S,I,...
    reservoirsID,...
    sensors,tau);

hsGSI_nom = GSI(h_nom_s,WLap,S,I,...
    reservoirsID,...
    sensors,tau);

disp(['RMSE(hreal,hGSI) = ' num2str(rmse(h,hGSI))]);

%% AW-GSI %%

[G_aw,B_aw,I2] = compute_GandB(graphInfo,hsGSI_nom);
[WLap2,Omega,Phi] = compute_LapfromGandB(hsGSI_nom,G_aw,B_aw,2,graphInfo);
Phi1 = diag(1./diag(Phi));

[fnomAW,~] = GSI(h_nom_s,WLap2,S,-I2,...
        graphInfo.reservoirsID,...
        sensors,tau);

res_AW = GSI_res(h_s-h_nom_s,WLap2,S,...
        graphInfo.reservoirsID,...
        sensors);
hAW = res_AW + fnomAW;

disp(['RMSE(hreal,hAW) = ' num2str(rmse(h,hAW))]);

% If running a nominal scenario, we change the used variables to their
% nominal analogues

if strcmp(scenario_type,'nominal') 
    h = h_nom;
    hGSI = hGSI_nom;
    hAW = fnomAW;
    demand = demand_nom;
    q = q_nom;
    demand_s = demand(idk);
    h_s = fix(h(sensors)*precision)/precision;
end

%% UKF configuration %%

initialStateGuess = hAW; % hAW is the default initial guess, using the AW-GSI state estimation

save temp_UKF_Measurement.mat sensors idk E invT 
save temp_UKF_State.mat alfa Phi1 Omega

yMeas = [h_s;demand_s];

ukf = unscentedKalmanFilter(...
    @StateFcn,... 
    @MeasurementFcn,... 
    initialStateGuess);

R = diag([0.0001*ones(length(sensors),1);0.0001*ones(length(idk),1)]);
ukf.MeasurementNoise = R;

ukf.ProcessNoise = 1*eye(length(hGSI));

Ku = 5; % Number of time steps in the UKF inner loop
xCorrectedUKF = zeros(Ku,N); % Corrected state estimates
PCorrected = zeros(Ku,N,N); % Corrected state estimation error covariances
rmse_h(:,1) = rmse(initialStateGuess,h);

%% UKF algorithm %%

K = 10;

sk=1;
for iter = 1:K
    for k=1:Ku
        disp(['### Iteration ' num2str(sk) ' out of ' num2str(Ku*K) ' ###'])
        [xCorrectedUKF(k,:), PCorrected(k,:,:)] = correct(ukf,yMeas);
        rmse_h(:,sk+1) = rmse(xCorrectedUKF(k,:),h);
        [~,~,JCorrectedUKF] = compute_qandB(xCorrectedUKF(k,:),E,T);

        predict(ukf);
        sk=sk+1;
    end
    %%% Comment the lines below to perform UKF-GSI (without state update) %%% 
    [Gaux,B,~] = compute_GandB(gI,xCorrectedUKF(end,:)');
    [WLap2,Omega,Phi] = compute_LapfromGandB(xCorrectedUKF(end,:)',Gaux,B,2,gI);
    Phi1 = diag(1./diag(Phi));
    save temp_UKF_State.mat alfa Phi1 Omega
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Plotting %

if ~ishandle(time)
    figure(time);
    set(gca, 'XLim', [0, get(gca, 'XLim') * [0; Ku*iter]])
    hold on
    yline(rmse(initialStateGuess,h),'k--','LineWidth',1.5)
end
figure(time);
hold on
fig = plot(0:length(rmse_h)-1,rmse_h,'.--','LineWidth',1.5,'MarkerSize',10);
ylabel('$RMSE(h^{est}_k,h^{real})$','interpreter','latex','fontsize',14)
xlabel('$k [\# iteration]$','interpreter','latex','fontsize',14); 
title('Head estimation error','interpreter','latex','fontsize',20); 
legend('AW-GSI baseline','UKF-AW-GSI','interpreter','latex','fontsize',12);

folder_name = ['test_leak' num2str(leak) '_' scenario_type];
folder = [erase(path,'\code') '/images/' folder_name '/'];
if ~exist(folder, 'dir')
    mkdir(folder)
end

saveas(fig,[folder 'time' num2str(time) '.fig'])
exportgraphics(gca, [folder 'time' num2str(time) '.png'], 'Resolution', 300)
close all
delete temp_UKF_State.mat temp_UKF_Measurement.mat

%% Functions %%

function x = StateFcn(x)
    load temp_UKF_State.mat alfa Omega Phi1
    x = alfa*x + (1-alfa)*Phi1*Omega*x;
end

function yk = MeasurementFcn(xk)
    load temp_UKF_Measurement.mat sensors idk E invT 
    yk(1:length(sensors)) = xk(sensors);
    [~,~,J] = compute_qandB(xk,E,diag(1./diag(invT)));
    yk(length(sensors)+1:length(sensors)+length(idk)) = -J(idk,:)*(invT*J'*xk).^(1/1.852);
end