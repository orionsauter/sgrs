%% PSSL - ILaRIS CVM
% Contains input parameters and initial conditions of the function
% "TMsys_Ball6_6DOF_v4.m" that describes the Actuator, Grabbing Finger (GF), Release Tip (RT), and Test Mass
% (TM) system according to paper "Prediction of the LISA-Pathfinder release
% mechanism in-flight performance" and the spring-tip concept. Integrates ODEs with the use of ODE45.

function TMrel_HG(varargin)
tic

%% Settings
if isempty(varargin)
    config_file = "config.xml";
else
    config_file = varargin{1};
end
config = readstruct(config_file);
fields = fieldnames(config);
for i = 1:numel(fields)
    % Assign values to variables. Using eval here is insecure, but where
    % would you even get a malicious config?
    eval(fields{i}+"="+config.(fields{i})+";")
end

timestamp = (num2str(fix(clock)));
savingid = timestamp((~isspace(timestamp)));
if exist("run_name", "var")
    savepath = run_name + "_" + savingid;
else
    savepath = savingid;
end

%% Parameters

% Retraction Speed
switch v_input
    case 'con'
        vrsx = 0;
        vrsy = 2e-6; % [m/s]
    case 'ts'
        load('v_input_2ums_array.mat') % retraction speed time-series
end

% Spring, Damping, & Friction
b=ksi*2*sqrt(k*mT);%0.01*2*sqrt(k*mT); % [N s/m] Damping constant

% Electrostatic stiffness
kGRS=-3.77e-9*Vinj^2-4.1e-9; % [N/m] from Giacomo's thesis

tsep1  = -timp; % [s] initial value for separation times
tsep2  = -timp; % [s]
tsept  = -timp; % [s]
imp = 0; % CGF impulse has not occured
imn = 0; % PGF impulse has not occured
imy = 0;

% Adhesion force
bond1 = 1; % starting value for TM to GF1 bond
bond2 = 1; % starting value for TM to GF2 bond
[ X1g1, X2g1, X3g1 ] = f_adh_predef( Ad_G1 );
[ X1g2, X2g2, X3g2 ] = f_adh_predef( Ad_G2 );
[ X1t1, X2t1, X3t1 ] = f_adh_predef( Ad_T1 );
[ X1t2, X2t2, X3t2 ] = f_adh_predef( Ad_T2 );

%% Initial Conditions

d = -(-xG10-0.5/k)/(1+(k/kT));%3.8e-3; % [m] Offset between GF and unstretched spring length

xT10=-k*d/kT;%-k*d/kT+sTM/2;%0; % [m] RT1 position
vT10=0; % [m/s] RT1 speed
xT20=k*d/kT;%k*d/kT-sTM/2;%0; % [m] RT2 position
vT20=0; % [m/s] RT2 speed 

y0 = [xG10 ; xG20 ; xT10 ; vT10 ; xT20 ; vT20 ; xTM0' ; vTM0' ; bTM0' ; wTM0'; bond1; bond2; imp; imn; imy; tsep1; tsep2; tsept];
jpatt = zeros([numel(y0),numel(y0)]);
jpatt(3,[4,11]) = 1;
jpatt(5,[6,11]) = 1;
jpatt(7,10) = 1;
jpatt(8,11) = 1;
jpatt(9,12) = 1;
jpatt(10,[3:9,13:15]) = 1;
jpatt(11,:) = 1;
jpatt(12,[3:9,13:15]) = 1;
jpatt(13,16) = 1;
jpatt(14,17) = 1;
jpatt(15,18) = 1;
jpatt(16,[3:9,13:15]) = 1;
jpatt(17,[3:9,13:15]) = 1;
jpatt(18,[3:9,13:15]) = 1;

%% Integration
% options = odeset("JPattern",jpatt,"Vectorized","on","Events",@unbond);
options = odeset("JPattern",jpatt,"Events",@(t,y)unbond(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG));
%[t,y] = ode45(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0);%, options);
ie = [];
t = [tspan(1), tspan(1)];
y = transpose(y0);
t_out = [];
y_out = [];
chkpts = 0:0.05:1.0;
chki = 1;
dur = tspan(end) - tspan(1);
t_events = [];
while t(end) < tspan(end)
    if tspan(1)/dur > chkpts(chki)
        fprintf("%.2f: %.0f%%\n",toc,chkpts(chki)*100);
        chki = chki + 1;
    end
    if numel(ie) > 0
        % Event occurred
        tspan(1) = min(te);
        t_events = cat(1, t_events, tspan(1));
        ye = ye(te<=tspan(1),:);
        te = te(te<=tspan(1),:);
        y0 = ye(1,:);
        [dydt, dLg1, Fg1, dLg2, Fg2, dLt1, Ft1, dLt2, Ft2] = TMsys_HG(tspan(1),y0',...
            ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,...
            vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,...
            I,T,Fx,Fz,tr,timp,rgf,muT,muG);
    else
        if t(1) > tspan(1) && t(end) < tspan(end)
            tspan(1) = t(end);
            y0 = y(1,:);
            fprintf("t = %f: Integration stopped\n",t(end));
            break
        else
            % Start of integration
            y0 = y(1,:);
        end
    end
    if any(ie == 1)
        e1 = find(ie==1,1);
        y0(:,21) = 1;
        y0(:,24) = te(e1);
        fprintf("%f: GF1/TM contact\n",te(e1));
    end
    if any(ie == 2)
        e2 = find(ie==2,1);
        y0(:,22) = 1;
        y0(:,25) = te(e2);
        fprintf("%f: GF2/TM contact\n",te(e2));
    end
    if any(ie == 3)
        e3 = find(ie==3,1);
        y0(:,19) = ~ye(e3,19);
        if ye(e3,21)==0
            y0(:,21) = 1;
            y0(:,24) = te(e3);
            fprintf("%f: GF1 bond change\n",te(e3));
        end
    end
    if any(ie == 4)
        e4 = find(ie==4,1);
        y0(:,20) = ~ye(e4,20);
        if ye(e4,22)==0
            y0(:,22) = 1;
            y0(:,25) = te(e4);
            fprintf("%f: GF2 bond change\n",te(e4));
        end
    end
    % RT to TM contact – Set RT vel to TM y
    if any(ie == 5)
        e5 = find(ie==5,1);
        y0(:,3) = y0(:,3) + dLt1;
        y0(:,4) = y0(:,11);
        fprintf("%f: RT1/TM contact\n",te(e5));
    end
    % RT to TM contact – Set RT vel to TM y
    if any(ie == 6)
        e6 = find(ie==6,1);
        y0(:,5) = y0(:,5) - dLt2;
        y0(:,6) = y0(:,11);
        fprintf("%f: RT2/TM contact\n",te(e6));
    end
    % RT to GF contact – Set RT vel to GF, and pos to full length
    if any(ie == 7)
        e7 = find(ie==7,1);
        y0(:,3) = ye(e7,1) - Lt;
        y0(:,4) = y0(:,1);
        fprintf("%f: GF1/RT1 contact\n",te(e7));
    end
    if any(ie == 8)
        e8 = find(ie==8,1);
        y0(:,5) = ye(e8,2) + Lt;
        y0(:,6) = y0(:,2);
        fprintf("%f: GF2/RT2 contact\n",te(e8));
    end
    y0(:,23) = y0(:,21) & y0(:,22);
    y0(:,26) = max(y0(:,24),y0(:,25)).*y0(:,23) + -timp.*~y0(:,23);
    [t,y,te,ye,ie] = ode45(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    % [t,y,te,ye,ie] = ode23s(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    % [t,y,te,ye,ie] = ode15s(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    if numel(ie) > 0
        % fprintf("%f: %d", te(1), ie(1));
        t_out = cat(1, t_out(1:end-1,:), t(t<=te(1),:));
        y_out = cat(1, y_out(1:end-1,:), y(t<=te(1),:));
    else
        t_out = cat(1, t_out(1:end-1,:), t);
        y_out = cat(1, y_out(1:end-1,:), y);
    end
    % [t,y,te,ye,ie] = ode15s(@(t,y)TMsys_HGvec(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
end

%% Output
% state vector y = [xG1 xG2 xT1 vT1 xT2 vT2 xTM yTM zTM v_xTM v_yTM v_zTM an_xTM an_yTM an_zTM om_xTM om_yTM om_zTM]
t = t_out;
y = y_out;
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

bond1 = y(:,19);
bond2 = y(:,20);
imp = y(:,21);
imn = y(:,22);
imy = y(:,23);
tsep1 = y(:,24);
tsep2 = y(:,25);
tsept = y(:,26);

clear y % saves memory

dLg1=xG1-xTM(:,2); % yG1-yTM
dLg2=xTM(:,2)-xG2; % yTM-yG2
dLt1=xT1-xTM(:,2); % yT1-yTM
dLt2=xTM(:,2)-xT2; % yTM-yT2


Fg1=kG*dLg1.*(dLg1<=0) + X1g1.*dLg1.*exp(-X2g1.*dLg1.^X3g1).*(dLg1>0 & bond1==1);
Fg2=kG*dLg2.*(dLg2<=0) + X1g2.*dLg2.*exp(-X2g2.*dLg2.^X3g2).*(dLg2>0 & bond2==1);
Ft1=kT*dLt1.*(dLt1<=0) + X1t1.*dLt1.*exp(-X2t1.*dLt1.^X3t1).*(dLt1>0); 
Ft2=kT*dLt2.*(dLt2<=0) + X1t2.*dLt2.*exp(-X2t2.*dLt2.^X3t2).*(dLt2>0);

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
    mkdir(savepath)
    copyfile(config_file,savepath);
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
L = legend('RT1', 'RT2', 'GF1', 'GF2', 'TM');
L.AutoUpdate = 'off';
xline(t_events);
title 'CM & TM Y-axis Displacement'
xlabel 'Time [sec]'
ylabel 'Displacement [\mum]'
grid on

figure(2) % TM Speed
plot(t,vTM/1e-6,'LineWidth',2)
L = legend('vTMx', 'vTMy', 'vTMz');
L.AutoUpdate = 'off';
xline(t_events);
title 'TM Speed'
xlabel 'Time [sec]'
ylabel 'Speed [\mum/s]'
grid on

figure(3) % TM Orientation
plot(t,bTM,'LineWidth',2)
L = legend('\psi', '\theta', '\phi');
L.AutoUpdate = 'off';
xline(t_events);
title 'TM Orientation'
xlabel 'Time [sec]'
ylabel 'Angle [rad]'
grid on

figure(4) % TM Ang. Speed
plot(t,wTM,'LineWidth',2)
L = legend('\omega_{psi}', '\omega_{theta}', '\omega_{phi}');
L.AutoUpdate = 'off';
xline(t_events);
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
L = legend('RT1', 'RT2');
L.AutoUpdate = 'off';
xline(t_events);
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
legend x z
xlabel 'Time [s]'
ylabel 'Displacement [\mum]'
grid on

if plots==1 && saveplot==1
    if savedata~=1
        cd(savepath)
    end
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
% cd 'C:\Users\anthony.davila\OneDrive - University of Florida\Shared Documents\Design\CVM\MATLAB Scripts'
% to = 'anthony.davila@ufl.edu';
% subject = 'Sim done';
% body = 'Yes';
% 
% sendolmail(to,subject,body)
runtime = toc
end

function [value,isterminal,direction] = unbond(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG)
    [dydt, dLg1, Fg1, dLg2, Fg2, dLt1, Ft1, dLt2, Ft2] = TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG);
    bond1 = y(19,1,:);
    bond2 = y(20,1,:);
    % Bonds are broken when exceeding the elongation threshold. To re-form a bond, contact must occur (dLg<=0).
    bond1 = 1*(dLg1<=0)+1*(dLg1>0 & dLg1<tr & bond1==1);
    bond2 = 1*(dLg2<=0)+1*(dLg2>0 & dLg2<tr & bond2==1);
    c1 = (y(1,1,:)-y(3,1,:)-Lt); % GF RT contact
    c2 = (y(5,1,:)-y(2,1,:)-Lt);
    value = [dLg1-tr; dLg2-tr; bond1; bond2; dLt1; dLt2; c1; c2];
    % isterminal = [1; 1; 1; 1; 0; 0; 0; 0];
    isterminal = [1; 1; 1; 1; 1; 1; 1; 1];
    direction = [1; 1; 0; 0; 0; 0; 0; 0];
end