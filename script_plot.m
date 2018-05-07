% Simple script to create a GR-FMD data fit and plot.
clear;

% Load in variables and define some constants.
load('Mag_data.mat');
Mc=1.3; % See scrip_Mc.m for the determination of 1.3 as the Mc.
dM=0.05;
GREY=[0.85,0.85,0.85];

% GR-FMD fit parameters from function call.
[b, b_err, a, R2_b,~,Mgr,Ngr,ngr]=Bval(M, Mc,dM);

% Define some values for plotting.
po=[-b,a];
Mgr_fit=Mc:dM:max(M);
Ngr_fit=10.^polyval(po,Mgr_fit);

% Output fit statistics to user prompt.
fprintf('b-value fit: %0.2f ± %0.2f  (R² %0.3f)\n',b,b_err, R2_b);

% Plot GR-FMD and fit.
figure(1); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlabel('Magnitude'); ylabel('Count');
xlim([0 5.5]); ylim([0.7 1.5*max(Ngr)]);