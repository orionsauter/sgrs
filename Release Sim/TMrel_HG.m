%% PSSL - ILaRIS CVM
% Contains input parameters and initial conditions of the function
% "TMsys_Ball6_6DOF_v4.m" that describes the Actuator, Grabbing Finger (GF), Release Tip (RT), and Test Mass
% (TM) system according to paper "Prediction of the LISA-Pathfinder release
% mechanism in-flight performance" and the spring-tip concept. Integrates ODEs with the use of ODE45.

clc, clear all, close all

%% Settings
global imp imn imy bond1 bond2 tsep1 tsep2 tsept

plots=1;
savedata=0;
saveplot=0;

tspan = [0 20]; % Time array [s]

Ad_G1=0; % GF1 adhesion force profile
Ad_G2=2; % GF2
Ad_T1=0; % RT1
Ad_T2=0; % RT2

v_input = 'con'; % 'con' for constant speed (faster simulation), 'ts' for time series speed profile (slower sim)

timestamp = (num2str(fix(clock)));
savingid = timestamp((~isspace(timestamp)));
savepath = ['C:\Users\anthony.davila\OneDrive - University of Florida\Shared Documents\Design\CVM\MATLAB Scripts\',savingid,'\'];

%% Parameters

% Time
ts1 = 50e-6; % [s] GF1 actuation delay time
ts2 = 0; % [s] GF2 actuation delay time
tp = 1; % [s] pause retraction of GFs
tc = 6; % [s] continue retraction of GFs

% Retraction Speed
switch v_input
    case 'con'
        vrsx = 0;
        vrsy = 2e-6; % [m/s]
    case 'ts'
        load('v_input_2ums_array.mat') % retraction speed time-series
end

% Dimensions
sTM = 30e-3; % [m] TM cube side length
sEH = 32e-3; % [m] EH cube side length
%d = 3.8e-3; % [m] Offset between GF and unstretched spring length
Lt = 0.5e-4; % [m] Maximum length that the RT extends from the GF
alp = 45*pi/180; % [rad] TM indent inclination
ac = 14.5e-3; %[m] distance of TM-GF contact point along y
bc = 2.5e-3; %[m] distance of TM-GF contact point ialong x or z
rgf = 4e-3/2; % [m] radius of GF tip

% Mass & Inertia
mG = (143/1000^3)*4430; % [kg] GF mass - annealed, grade 5 Titanium
mT = (4*pi*2.5^2/1000^3)*3200; % [kg] RT mass - Silicon Nitride
M = 0.5378; % [kg] ILaRIS TM mass - 70Au-30Pt
I = eye(3)*(1/6)*M*sTM^2; % [kg-m^2] TM matrix of inertia

% Spring, Damping, & Friction
k=1.317e4;%0.75*(4.45/0.0254);%1.317e4; % [N/m] Spring constant
kG=5e5;%6e7;%114e9*8.6/(16*1000); % [N/m] GF equivalent spring constant (from GF prelim geometry)
kT=5e5;%1e7;%1.515e7; %(pi*2.5e-3^2)*290e9/4e-3 % [N/m] RT equivalent spring constant
ksi=0.8; % Damping ratio
b=ksi*2*sqrt(k*mT);%0.01*2*sqrt(k*mT); % [N s/m] Damping constant
muT = 0.1; % Static friction coefficient of silicon nitride ceramic tips
muG = 0.1 ; % Static friction coefficient of 70Au-30Pt GFs (assumed)

% Electrostatic stiffness
Vinj=4; %[V] Injection signal bias
kGRS=-3.77e-9*Vinj^2-4.1e-9; % [N/m] from Giacomo's thesis

% Impulses (from GPRM data)
ixp = 75e-6; % [kg m/s] TM release x impulse on CGF side
ixn = 85e-6; % [kg m/s] TM release x impulse on PGF side
izp = -40e-6; % [kg m/s] TM release z impulse on CGF side
izn = -20e-6; % [kg m/s] TM release z impulse on PGF side
iy = 60e-6; % [kg m/s] TM release residual y impulse
Fx = 0; % [N] Lateral force x
Fz = 0;% [N] Lateral force y
T = [0,0,0]; % [N-m] Pure torque
timp = 0.5; % [s] time interval in which to apply lateral impulse

tsep1  = -timp; % [s] initial value for separation times
tsep2  = -timp; % [s]
tsept  = -timp; % [s]
imp = 0; % CGF impulse has not occured
imn = 0; % PGF impulse has not occured
imy = 0;

% Adhesion force
tr = 4e-6; % [m] bond elongation threshold (if exceeded the adhesion bond is broken)
bond1 = 1; % starting value for TM to GF1 bond
bond2 = 1; % starting value for TM to GF2 bond
[ X1g1, X2g1, X3g1 ] = f_adh_predef( Ad_G1 );
[ X1g2, X2g2, X3g2 ] = f_adh_predef( Ad_G2 );
[ X1t1, X2t1, X3t1 ] = f_adh_predef( Ad_T1 );
[ X1t2, X2t2, X3t2 ] = f_adh_predef( Ad_T2 );

%% Initial Conditions

xG10=-3.26e-6;%-3.26e-6+sTM/2;%-200/kG;% % [m] GF1 position
vG10=0; % [m/s] GF1 speed
xG20=3.26e-6;%3.26e-6-sTM/2;%200/kG; % [m] GF2 position
vG20=0; % [m/s] GF1 speed

d = -(-xG10-0.5/k)/(1+(k/kT));%3.8e-3; % [m] Offset between GF and unstretched spring length

xT10=-k*d/kT;%-k*d/kT+sTM/2;%0; % [m] RT1 position
vT10=0; % [m/s] RT1 speed
xT20=k*d/kT;%k*d/kT-sTM/2;%0; % [m] RT2 position
vT20=0; % [m/s] RT2 speed 

xTM0=[0,0,0]; % [m] TM position 
vTM0=[0,0,0]; % [m/s] TM velocity
bTM0=[0,0,0]; % [rad] TM orientation (body system)
wTM0=[0,0,0]; % [rad/s] TM angular velocity (body system)

y0 = [xG10 ; xG20 ; xT10 ; vT10 ; xT20 ; vT20 ; xTM0' ; vTM0' ; bTM0' ; wTM0'];

%% Integration

% options = odeset('AbsTol',1e-2);
[t,y] = ode45(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0);%, options);

%% Output
% state vector y = [xG1 xG2 xT1 vT1 xT2 vT2 xTM yTM zTM v_xTM v_yTM v_zTM an_xTM an_yTM an_zTM om_xTM om_yTM om_zTM]

xG1=y(:,1);
xG2=y(:,2);
xT1=y(:,3);
vT1=y(:,4);
xT2=y(:,5);
vT2=y(:,6);
xTM=y(:,7:9);
vTM=y(:,10:12);
bTM=y(:,13:15);
wTM=y(:,16:18);

clear y % saves memory

dLg1=xG1-xTM(:,2); % yG1-yTM
dLg2=xTM(:,2)-xG2; % yTM-yG2
dLt1=xT1-xTM(:,2); % yT1-yTM
dLt2=xTM(:,2)-xT2; % yTM-yT2

Fg1 = zeros(size(t));
Fg2 = zeros(size(t));
Ft1 = zeros(size(t));
Ft2 = zeros(size(t));
for ii = 1:size(t,1)
bond1 = 1*(dLg1(ii)<=0)+1*(dLg1(ii)>0 && dLg1(ii)<tr && bond1==1);
bond2 = 1*(dLg2(ii)<=0)+1*(dLg2(ii)>0 && dLg2(ii)<tr && bond2==1);
Fg1(ii)=kG*dLg1(ii)*(dLg1(ii)<=0) + X1g1*dLg1(ii).*exp(-X2g1*dLg1(ii).^X3g1)*(dLg1(ii)>0 && bond1==1);
Fg2(ii)=kG*dLg2(ii)*(dLg2(ii)<=0) + X1g2*dLg2(ii).*exp(-X2g2*dLg2(ii).^X3g2)*(dLg2(ii)>0 && bond2==1);
Ft1(ii)=kT*dLt1(ii)*(dLt1(ii)<=0) + X1t1*dLt1(ii).*exp(-X2t1*dLt1(ii).^X3t1)*(dLt1(ii)>0); 
Ft2(ii)=kT*dLt2(ii)*(dLt2(ii)<=0) + X1t2*dLt2(ii).*exp(-X2t2*dLt2(ii).^X3t2)*(dLt2(ii)>0);
end

%% GF Speed Approximate (to confirm input is well applied)

vG1e=diff(xG1)./diff(t);
vG2e=diff(xG2)./diff(t);

%% TM Corner Positions

psi=bTM(end,1);
the=bTM(end,2);
phi=bTM(end,3);

rEH=sEH/2*[-1,-1,-1 ; 1,-1,-1 ; 1,-1,1 ; -1,-1,1 ; -1,1,-1 ; 1,1,-1 ; 1,1,1 ; -1,1,1]; % [m] EH corners
rc=sTM/2*[-1,-1,-1 ; 1,-1,-1 ; 1,-1,1 ; -1,-1,1 ; -1,1,-1 ; 1,1,-1 ; 1,1,1 ; -1,1,1]; % [m] corner position from body system
R=[cos(the)*cos(phi) , cos(the)*sin(phi) , -sin(the) ; sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi) , sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) , sin(psi)*cos(the) ; cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi) , cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) , cos(psi)*cos(the)];
rcf=(R*rc')'; % [m] corner position from fixed system
rcfi=rcf+ones(size(rcf)).*(xTM(end,:)); % [m] TM corners final position

%% Contact conditions (RT correction)

cG_RT1=xG1-xT1>=Lt;
cG_RT2=xT2-xG2>=Lt;
cG_TM1=abs(xG1-xTM(:,2))<tr;
cG_TM2=abs(xG2-xTM(:,2))<tr;
cRT_TM1=abs(xT1-xTM(:,2))<tr;
cRT_TM2=abs(xT2-xTM(:,2))<tr;

%% Output messages

fprintf('\n\nOutput\n\n Max. TM Residual Velocity = %.3d m/s\n',max(abs(vTM(:,2))))
fprintf(' Final TM Residual Velocity = %.3d m/s\n',abs(vTM(end,2)))
if abs(vTM(end,2))>=4.5e-6     % Max residual velocity from "ResVel_2DOF.m"
    fprintf(' TM collision with EH\n')
else
    fprintf(' TM stopped by actuation with SF = %.1f\n',4.5e-6/abs(vTM(end,2)))
end

% Saved report

if savedata==1 || saveplot==1
        mkdir(savingid)
end

if savedata==1
    cd(savepath)
    save(['TMrel_',savingid,'.mat'],'t','xG1','xG2','xT1','vT1','xT2','vT2','xTM','vTM','bTM','wTM','-v7.3')
    fileID = fopen(['Output_',savingid,'.txt'],'w');
    fprintf(fileID,'Input\n\n Adhesion Set:\n  GF1 = %.0f\n  GF2 = %.0f\n  RT1 = %.0f\n  RT2 = %.0f \n',Ad_G1,Ad_G2,Ad_T1,Ad_T2);
    fprintf(fileID,' Retraction start:\n  GF1 = %.2d s\n  GF2 = %.2d s \n',ts1,ts2);
    fprintf(fileID,' System parameters:\n  k = %.2d m/s\n  kGF = %.2d m/s\n  kRT = %.2d m/s\n  kGRS = %.2d m/s\n  ksi = %.2d\n  b = %.2d N-s/m\n',k,kG,kT,kGRS,ksi,b);
    fprintf(fileID,'  spring_offset (d) = %.2d m\n  RT_to_GF_limit (Lt) = %.2d m\n',d,Lt);
    fprintf(fileID,' Initial conditions:\n  GF1: x0 = %.2d m , v0 = %.2d m/s\n  GF2: x0 = %.2d m , v0 = %.2d m/s\n  RT1: x0 = %.2d m , v0 = %.2d m/s\n  RT2: x0 = %.2d m , v0 = %.2d m/s\n  TM: x0 = %.2d m , y0 = %.2d m , z0 = %.2d m , vx0 = %.2d m/s , vy0 = %.2d m/s , vz0 = %.2d m/s\n',xG10,vG10,xG20,vG20,xT10,vT10,xT20,vT20,xTM0,vTM0);
    fprintf(fileID,'\n\nOutput\n\n Max. TM Residual Velocity = %.3d m/s\n',max(abs(vTM(:,2))));
    fprintf(fileID,' Final TM Residual Velocity = %.3d m/s\n',abs(vTM(end,2)));
    if abs(vTM(end,2))>=4.5e-6     % Max residual velocity from "ResVel_2DOF.m"
        fprintf(fileID,' TM collision with EH');
    else
        fprintf(fileID,' TM stopped by actuation with SF = %.1f\n',4.5e-6/abs(vTM(end,2)));
    end
    fclose(fileID);
end

%% Plots
if plots==1

% figure(1) % Position
% mins = min(abs([xT1,xT2,xG1,xG2])); % obtains minimum offset to TM center (for plot visualization)
% plot(t,(xT1-mins(1))/1e-6,t,(xT2+mins(2))/1e-6,t,(xG1-mins(3))/1e-6,t,(xG2+mins(4))/1e-6,t,xTM(:,2)/1e-6,'LineWidth',2)
% legend RT1 RT2 GF1 GF2 TM
% title 'CM & TM Y-axis Displacement'
% xlabel 'Time [sec]'
% ylabel 'Displacement [\mum]'
% grid on

figure(1) % Position
plot(t,xT1/1e-6,t,xT2/1e-6,t,xG1/1e-6,t,xG2/1e-6,t,xTM(:,2)/1e-6,'LineWidth',2)
legend RT1 RT2 GF1 GF2 TM
title 'CM & TM Y-axis Displacement'
xlabel 'Time [sec]'
ylabel 'Displacement [\mum]'
grid on

figure(2) % TM Speed
plot(t,vTM/1e-6,'LineWidth',2)
legend vTMx vTMy vTMz
title 'TM Speed'
xlabel 'Time [sec]'
ylabel 'Speed [\mum/s]'
grid on

figure(3) % TM Orientation
plot(t,bTM,'LineWidth',2)
legend \psi \theta \phi
title 'TM Orientation'
xlabel 'Time [sec]'
ylabel 'Angle [rad]'
grid on

figure(4) % TM Ang. Speed
plot(t,wTM,'LineWidth',2)
legend \omega_{psi} \omega_{theta} \omega_{phi}
title 'TM Angular Speed'
xlabel 'Time [sec]'
ylabel 'Angular Speed [rad/s]'
grid on

figure(5) % GF forces
plot(t,Fg1,t,Fg2,'LineWidth',2)
legend GF1 GF2
title 'Force between GFs and TM'
xlabel 'Time [sec]'
ylabel 'Force [N]'
grid on

figure(6) % RT forces
plot(t,Ft1,t,Ft2,'LineWidth',2)
legend RT1 RT2
title 'Force between RTs and TM'
xlabel 'Time [sec]'
ylabel 'Force [N]'
grid on

figure(7) % Contact
    subplot(3,1,1)
    plot(t,cG_RT1,t,cG_RT2,'LineWidth',1)
    xlim([t(1)-0.1 t(end)+0.1])
    ylim([-0.1 1.1])
    legend GF-RT_1 GF-RT_2
    title('GF to RT contact')

    subplot(3,1,2)
    plot(t,cG_TM1,t,cG_TM2,'LineWidth',1)
    xlim([t(1)-0.1 t(end)+0.1])
    ylim([-0.1 1.1])
    legend GF-TM_1 GF-TM_2
    title('GF to TM contact')
    ylabel 'Contact (1 = Yes, 0 = No)'

    subplot(3,1,3)
    plot(t,cRT_TM1,t,cRT_TM2,'LineWidth',1)
    xlim([t(1)-0.1 t(end)+0.1])
    ylim([-0.1 1.1])
    legend RT-TM_1 RT-TM_2
    title('RT to TM contact')
    xlabel 'Time [sec]'
    
figure(8) % GF Speed Approximate
plot(t(1:numel(vG1e)),vG1e,t(1:numel(vG2e)),-vG2e,'linewidth',2);
title 'Grabbing Fingers Input Speed'
xlabel 'Time [s]'
ylabel 'Speed [m/s]'
legend vGF1 -vGF2
grid on
end

figure(9) % 3D TM position in EH
idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]'; % cube side index
xc = rEH(:,1); yc = rEH(:,2); zc = rEH(:,3);
xc1 = rc(:,1); yc1 = rc(:,2); zc1 = rc(:,3);
xc2 = rcfi(:,1); yc2 = rcfi(:,2); zc2 = rcfi(:,3);
patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.1);
hold on
patch(xc1(idx), yc1(idx), zc1(idx), 'b', 'facealpha', 0.1);
hold on
patch(xc2(idx), yc2(idx), zc2(idx), 'g', 'facealpha', 0.1);
view(3);
legend EH TM_i TM_f
title 'TM positions with respect to the EH'
axis equal
xlabel 'X [m]'
ylabel 'Y [m]'
zlabel 'Z [m]'

figure(10) % TM position in X and Z axes
plot(t,xTM(:,1)/1e-6,t,xTM(:,3)/1e-6,'linewidth',2)
title 'TM position'
legend x y z
xlabel 'Time [s]'
ylabel 'Displacement [\mum]'
grid on

if plots==1 && saveplot==1
        cd(savepath)    
        saveas(figure(1),['x_',savingid,'.png'])
        saveas(figure(2),['vTM_',savingid,'.png'])
        saveas(figure(3),['bTM_',savingid,'.png'])
        saveas(figure(4),['wTM_',savingid,'.png'])
        saveas(figure(5),['Fgf_',savingid,'.png'])
        saveas(figure(6),['Frt_',savingid,'.png'])
        saveas(figure(7),['Cont_',savingid,'.png'])
        saveas(figure(8),['vGF_',savingid,'.png'])   
        saveas(figure(10),['TMpos_',savingid,'.png'])   
end

%% Send Email once it's done
cd 'C:\Users\anthony.davila\OneDrive - University of Florida\Shared Documents\Design\CVM\MATLAB Scripts'
to = 'anthony.davila@ufl.edu';
subject = 'Sim done';
body = 'Yes';

sendolmail(to,subject,body)