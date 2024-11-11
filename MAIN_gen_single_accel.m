clc; clearvars; close all;
set(groot,'defaultFigureColor','w')
set(groot,'defaultLineLineWidth',1)
g = 9.81; %m/sec2

% This code generates target spectrum-compatible fully-nonstationary 
% artificial seismic ground motions using the Continuous Wavelet Transform
% and a seed record.

% Further details are provided in the following document:
%
% Yanni H., Fragiadakis M., and Mitseas I.P. 
% "Wavelet-based stochastic model for the generation of fully
% non-stationary bidirectional seismic accelerograms". 
% Earthquake Engineering and Structural Dynamics. In review.

% Version 1.0 created by Hera Yanni, first release: 11th of November, 2024. 

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f_0 = lowest frequency (Hz) of the frequency range of the generated signals 
% fmax = highest frequency (Hz) of the frequency range of the generated signals
% rec_aG = the seed record in time and acceleration columns format 
% N = number of harmonics to be superposed
% nrecs = number of ground motions to be generated, 1
% N_iter = number of corrective iterations
% targSpec = type of target spectrum to be used

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the number of ground motions to be generated
nrecs = 1; % should be 1 for one record

% Load the seed record
rec_aG =  load('Kozani_1995_L.dat'); 
% Always convert the acceleration values to m/s^2 !!!
% e.g. Lefkada_1 and Lefkada_2 are in cm/sec^2
% Kozani has values in g
% Syntagma has values in m/s^2
ut_G1 = rec_aG(:,2).*g;     % or.*0.01 for cm/sec^2

% Define the desired frequency range of the ground motions
f_0 = 1/(2*pi); % lowest frequency (Hz), can't be less than 0.36/(2*pi)
fmax = 125/(2*pi); % highest frequency (Hz)

% Define the number of harmonics to be superposed
N = 500; 

% Define the desired number of corrective iterations
N_iter = 3;

% For targSpec = 'GMM', the model uses a GMM to generate the spectra and
% the accelerograms
% For targSpec = 'EC8', the model uses the EC8 spectrum to generate the 
% spectra and the accelerograms
% targSpec = 'GMM';
targSpec = 'EC8';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seed record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = rec_aG(:,1); % time
dt = t(3)-t(2); % time step
L = t(end);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the frequency/period range for target spectrum generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 1/dt;   % sampling frequency
Fn = Fs/2;   % Nyquist frequency

omega_max = 2*pi*Fn; % maximum frequency
omega_m = 2*pi*fmax; % rad/sec

omega_u = min(2*pi*fmax, omega_max); % rad/sec
omega_0 = 2*pi*f_0; % rad/sec

% Cut-off frequency upper limit check
if omega_m > omega_max
    msgbox('Maximum frequency must be less than the Nyquist frequency')
      return
end

% Cut-off frequency lower limit check
if omega_0 < 0.36 
    msgbox('Minimum frequency must be larger than 0.36 rad/s')
      return
end

omega = linspace(omega_0,omega_u,N)'; % frequencies for the harmonics
dom = mean(diff(omega)); % frequency domain integration step
T = flip(2*pi./omega); % periods for the target spectrum spectral values

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the target spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GMM
% The median and logarithmic standard deviation values are obtained from 
% the BSSA14 NGA-West2 model, based on the following:
%
% Boore DM, Stewart JP, Seyhan E, Atkinson GM. 
% NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for 
% Shallow Crustal Earthquakes. Earthquake Spectra. 2014;30:1057-1085.
% https://doi.org/10.1193/070113EQS184M.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EC8
% EC8 (EN1998-1-2004) code spectrum 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(targSpec,'EC8')
% %EC8 elastic spectrum
agR = 0.24; %g
gamma_I=1.00;
ag = agR*gamma_I*g; %m/s^2
soil = 'B'; % EC8 soil type
S=1.2; % S soil factor
PGA = ag*S/g;
zeta= 0.05;

% Target spectrum: EC8 spectrum
[Sa_mean,~,~] = EC8spectrumElastic(ag,soil,T,zeta); % m/s^2

elseif strcmp(targSpec,'GMM')

% GMM BSSA 14 
% Target response spectrum
Mw = 7; % Mw = Moment Magnitude
Rjb = 15; % Rjb = Joyner-Boore distance (km)
% Fault_Type    = 0 for unspecified fault
%               = 1 for strike-slip fault
%               = 2 for normal fault
%               = 3 for reverse fault
% normal fault
Fault_Type = 1; 
% Vs30 = shear wave velocity averaged over top 30 m in m/s
Vs30 = 800;
% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China or Turkey
%               = 4 for Italy
region = 4;
% z1            = Basin depth (km); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = 999 if unknown
z1 = 999;

% PGA
[PGA sigma_ag] = gmpe_bssa_2014(Mw, 0, Rjb, Fault_Type, region, z1, Vs30);

ag = PGA*g;
zeta = 0.05;

% Ground Motion Model
for i=1:N
    [Sa_mean(i) sigma_ln(i)] = gmpe_bssa_2014(Mw, T(i), Rjb, Fault_Type, region, z1, Vs30);
end

Sa_mean = Sa_mean.*g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user inputs. 
% Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Recorded Accelerogram response spectrum

[Sa_rec]=ARS([t'; ut_G1'],zeta,T');
Sa_rec = Sa_rec(4,:)';

%% Produce the spectrum compatible accelerogram

Ts = L; % accelerogram duration
nrecs=1;
accelerograms = zeros(nrecs,length(t));

% spectrum matching period range
T1=0.01:0.01:4;

if strcmp(targSpec,'EC8')
[Sa_,~,~] = EC8spectrumElastic(ag,soil,T1,zeta); % m/s^2

elseif strcmp(targSpec,'GMM')
for iij=1:length(T1)
     [Sa_(iij) ~] = gmpe_bssa_2014(Mw, T1(iij), Rjb, Fault_Type, region, z1, Vs30);
     Sa_(iij) = Sa_(iij)*g;
end
end

for i=1:nrecs
%% Stationary accelrogram    

% generate the stationary accelerogram
[a_G,G] = Stat_Accel_Gen(zeta,N,omega,Sa_mean,Ts,dom,t,Sa_,PGA,T1);

% Stationary accelerogram spectrum
[Sa_aG]=ARS([t';a_G'],zeta,T');
Sa_aG = Sa_aG(4,:)';

% CWT of the produced stationary artificial accelerogram
% Create filterbank
% Morlet wavelet = 'amor', Morse wavelet = 'morse'
f0=omega_0/(2*pi);
fu=omega_u/(2*pi);
fb_amor = cwtfilterbank_L2('SignalLength',numel(a_G),'SamplingFrequency',1/dt,'FrequencyLimits',[f0 fu],'VoicesPerOctave',48,'Wavelet','amor');
%fb_amor = cwtfilterbank_L2('SignalLength',numel(a_G),'SamplingFrequency',1/dt,'FrequencyLimits',[f0 fu],'VoicesPerOctave',48,'Wavelet','morse')

[wt_aGartif,f,~]=cwt_L2(a_G,'Filterbank',fb_amor);

[~,~,scls]=scales(fb_amor); % get the scales

%% Seed record

% CWT of the seed record
[wt_rec,f_rec,~]=cwt_L2(ut_G1,'Filterbank',fb_amor);

%% Produce the new accelerogram

% parameters of Eq. 17
max_pga_rec=max(abs(ut_G1),[],'all');
max_pga_artif=max(abs(a_G),[],'all');

lambda = max_pga_artif./max_pga_rec; % Eq. 17

% parameter of Eq. 18
mod_aG = abs(wt_aGartif);

% parameters of Eq. 16
mod_rec = abs(wt_rec);
w_max=max(mod_rec,[],'all'); 

wt_norm = (mod_rec./w_max); % Eq. 16

% Produce the new accelerogram's CWT coefficients
for ik=1:length(f)
    E_aG(ik) = trapz(t,mod_aG(ik,:)); % Eq.18
    E_rec(ik) = trapz(t,mod_rec(ik,:)); % Eq 19

    WT(ik,:) = (E_aG(ik)./(lambda*E_rec(ik))).*wt_norm(ik,:).*wt_aGartif(ik,:); % Eq. 20
end

% Inverse CWT 
ia_G = 0;
ia_G = icwt_L2(WT,scls,'amor',f,[f0+(dom/2*pi) fu]); % the new artificial accelerogram
%ia_G = icwt_L2(WT,scls,'morse',f,[f0+(dom/2*pi) fu]);

[PSa_iaG]=ARS([t'; ia_G;],zeta,T1);

[ut_G,PSa_aG]=FT_iter([t'; ia_G],N_iter,[T1;Sa_],T1,zeta,PGA*g);

% Baseline correction: 
% L: linear, Q: quadratic, C: cubic
[ut_G] = baseline_correction(ut_G','C');

accelerograms(i,:) = ut_G;
[PSa(i,:)]=PSa_aG(N_iter,:);

% close all
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write the results 
acc1=accelerograms(1,:);
 acc_ccl = [t, acc1']
 writematrix(acc_ccl,'acc_ccl.txt','Delimiter',' ') 

 %  save('a_G.mat','a_G');
 %  save('wt_aG.mat','wt_aG');
 %  save('f.mat','f');
 %  save('wt_rec.mat','wt_rec');
 % save('f_rec.mat','f_rec');
 %  save('WT.mat','WT');
 %  save('ia_G.mat','ia_G');
 %  save('wt_norm.mat','wt_norm'); 
 %  save('ut_G.mat','ut_G');
 %  save('wt_aGartif.mat','wt_aGartif');

 %% PLOTS

 % Plot all the information for the stationary accelerogram 
figure()
% PSD
subplot(2,2,1)
hold on; grid on; box on;
plot(G.*10^4,omega,'b','LineWidth',3) % PSD in cm2/sec3
ylim([0,omega_u])
xlim([0,max(G.*10^4)])
set(gca,'XDir','reverse','FontSize',21,'FontName','Times New Roman')
xlabel('$G(\omega)$ [g]','FontSize',24,'interpreter','latex');
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
% CWT
subplot(2,2,2)
hold on; grid on; box on;  
ylim([omega_0,omega_u])
xlim([0,L]);
pcolor(t,2*pi.*f,abs(wt_aGartif));
shading interp
colormap('jet')
set(gca,'FontSize',21,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
xlabel('Time [s]','FontSize',24,'interpreter','latex');
s.EdgeColor = 'none';
% spectrum vs target spectrum
subplot(2,2,3)
hold on; grid on; box on;
set(gca,'FontSize',21,'FontName','Times New Roman')
plot(T,Sa_mean./g,'r-.','LineWidth',3);
 plot(T,Sa_aG./g,'b','LineWidth',2);
xlabel('Period [s]','FontSize',24,'interpreter','latex');
ylabel('Spec. accel. [g]','FontSize',24,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}^{s}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','northeast')
xlim([0, 4])
%stationary artificial accelerogram
subplot(2,2,4)
title('Stationary accelerogram')
hold on; grid on; box on;
plot(t,a_G./g,'b')
set(gca,'FontSize',21,'FontName','Times New Roman')
yline(0,'k')
xlim([0, L])
ylim([-0.3, 0.301])
xlabel('Time [s]','FontSize',24,'interpreter','latex');
ylabel('$\alpha_s(t)$ [g]','FontSize',24,'interpreter','latex');

% Plot all the information for the seed record
figure()
subplot(2,2,1)
hold on; grid on; box on;
[fft1,FAm1] = FFTp(t(3)-t(2),ut_G1); 
plot(FAm1./g,2*pi*fft1,'b') 
ylim([0,omega_u])
xlim([0,0.1])
yticks([0 20 40 60 80 100 120])
yticklabels({'0','20','40','60','80','100','120'})
set(gca,'XDir','reverse','FontSize',21,'FontName','Times New Roman')
xlabel('$|F(\omega)|$ [g]','FontSize',24,'interpreter','latex');
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');

subplot(2,2,2)
hold on; grid on; box on;  
ylim([omega_0,omega_u])
xlim([0,L]);
pcolor(t,2*pi.*f_rec,abs(wt_rec));
shading interp
colormap('jet')
set(gca,'FontSize',21,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
xlabel('Time [s]','FontSize',24,'interpreter','latex');
s.EdgeColor = 'none';

subplot(2,2,3)
hold on; grid on; box on;
set(gca,'FontSize',21,'FontName','Times New Roman')
plot(T,Sa_mean./g,'r-.','LineWidth',3);
plot(T,Sa_rec./g,'b','LineWidth',2);
xlabel('Period [s]','FontSize',24,'interpreter','latex');
ylabel('Spec. accel. [g]','FontSize',24,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}^{rec}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','northeast')
xlim([0, 4])

subplot(2,2,4)
title('Seed record')
hold on; grid on; box on;
plot(t,ut_G1./g,'b')
set(gca,'FontSize',21,'FontName','Times New Roman')
yline(0,'k')
xlim([0, L])
ylim([-0.3, 0.301])
xlabel('Time [s]','FontSize',24,'interpreter','latex');
ylabel('$\alpha_r(t)$ [g]','FontSize',24,'interpreter','latex');

% plot the time-frequency modulating function
figure()
title('Eq. 16: time-frequency modulating function')
hold on; grid on; box on;
ylim([omega_0,omega_u])
xlim([0,L]);
zlim([0,max(abs(wt_norm),[],'all')]);
surf(t,2*pi.*f_rec,abs(wt_norm));
view(45,45)
shading interp
colormap('jet')
set(gca,'FontSize',19,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',22,'interpreter','latex');
xlabel('Time [s]','FontSize',22,'interpreter','latex');
zlabel('$\Phi(\omega,t)$','FontSize',22,'interpreter','latex');
s.EdgeColor = 'none';
hcb = colorbar('location','EastOutside');
% print('TF_Envelope.png','-dpng')

% PLOT THE PRODUCED ACCELEROGRAM
figure
subplot(3,3,[1,4])
hold on; grid on; box on;
[fft1,FAm1] = FFTp(t(3)-t(2),acc1); 
plot(FAm1./g,2*pi*fft1,'b') 
ylim([0,omega_u])
xlim([0,0.3])
set(gca,'XDir','reverse','FontSize',21,'FontName','Times New Roman')
xlabel('$|F(\omega)|$ [g]','FontSize',24,'interpreter','latex');
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');

subplot(3,3,[2,3,5,6])
hold on; grid on; box on;  
ylim([omega_0,omega_u])
xlim([0,L]);
[wt_f1,f1]=cwt_L2(acc1,'Filterbank',fb_amor);
pcolor(t,2*pi.*f1,abs(wt_f1));
shading interp
colormap('jet')
set(gca,'FontSize',21,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
xlabel('Time [s]','FontSize',24,'interpreter','latex');
s.EdgeColor = 'none';

subplot(3,3,7)
hold on; grid on; box on;
set(gca,'FontSize',21,'FontName','Times New Roman')
plot([0 T1],[PGA Sa_./g],'r-.','LineWidth',3);
plot(T1,PSa(1,:)./g,'b','LineWidth',2);
xlabel('Period [s]','FontSize',24,'interpreter','latex');
ylabel('Spec. accel. [g]','FontSize',24,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','northeast')
xlim([0, 4])
    
subplot(3,3,[8,9])
title('Produced accelerogram')
hold on; grid on; box on;
plot(t,acc1./g,'b')
set(gca,'FontSize',21,'FontName','Times New Roman')
yline(0,'k')
%xlim([0, L])
xlim([0, L])
ylim([-0.4, 0.401])
xlabel('Time [s]','FontSize',24,'interpreter','latex');
ylabel('$\alpha(t)$ [g]','FontSize',24,'interpreter','latex');
 yline(PGA,'r--','LineWidth',2);
 yline(-PGA,'r--','LineWidth',2);
 % printScript
