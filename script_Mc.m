% Simple script to create find Mc from GR-FMD data fits.  Based on the
% method of Wiemer & Wyss (2000). Bull. Seismol. Soc. Am. 90,859?869, 
% doi:10.1785/0119990114.
clear;

% Load in variables and define some constants.
load('Mag_data.mat');
Mc_list=0:0.01:3.2;
dM=0.05;
GREY=[0.85,0.85,0.85];

% Preallocate space in memory.
R=zeros(size(Mc_list));
X=R;
b=R;

% Loop through every Mc truncation value and fit GR-FMD.
for i=1:length(Mc_list)
    [temp,~,~,temp2, temp3, ~,~,~]=Bval(M, Mc_list(i), dM);
    b(i)=temp;
    R(i)=temp2;
    X(i)=temp3;
end;

% Some error handling.
R(R>1)=NaN;
R(R<0)=NaN;

% Plot GR-FMD fit statistics, visually determine fit stability and Mc.
figure(2); clf;
subplot(211); plot(Mc_list, log10(1-R)); xlabel('Magnitude'); ylabel('log_{10}(1-R^2)');
subplot(212); plot(Mc_list, b); xlabel('Magnitude'); ylabel('b-value');