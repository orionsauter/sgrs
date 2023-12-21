%% PSSL - ILaRIS CVM
% Contains input parameters and initial conditions of the function
% "TMsys_Ball6_6DOF_v4.m" that describes the Actuator, Grabbing Finger (GF), Release Tip (RT), and Test Mass
% (TM) system according to paper "Prediction of the LISA-Pathfinder release
% mechanism in-flight performance" and the spring-tip concept. Integrates ODEs with the use of ODE45.

function [captured] = TMrel_HG(varargin)
tic

%% Settings
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'config','config.xml');
addOptional(p,'estimateGF',0);
addOptional(p,'useEnergy',0);
parse(p,varargin{:});
config = readstruct(p.Results.config);
% Hard code the field names to get the order right
fields = {...
    'run_name'; 'plots'; 'savedata'; 'saveplot'; 'savestats'; 'tspan';...
    'Ad_G1'; 'Ad_G2'; 'Ad_T1'; 'Ad_T2'; 'v_input'; 'ts1'; 'ts2';...
    'tp'; 'tc'; 'sTM'; 'sEH'; 'd'; 'Lt'; 'alp'; 'ac'; 'bc'; 'rgf';...
    'mG'; 'mT'; 'M'; 'I'; 'k'; 'kG'; 'kT'; 'ksi'; 'muT'; 'muG'; 'Vinj';...
    'ixp'; 'ixn'; 'izp'; 'izn'; 'iy'; 'Fx'; 'Fz'; 'T'; 'timp'; 'tr';...
    'xG10'; 'vG10'; 'xG20'; 'vG20'; 'xTM0'; 'vTM0'; 'bTM0'; 'wTM0'; 'vrsy'};
for i = 1:numel(fields)
    % Assign values to variables. Using eval here is insecure, but where
    % would you even get a malicious config?
    if isfield(config, fields{i})
        eval(fields{i}+"="+config.(fields{i})+";");
    end
    % Arguments override config
    if isfield(p.Unmatched, fields{i})
        eval(fields{i}+"="+p.Unmatched.(fields{i})+";");
    end
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
        if ~exist('vrsy','var')
            vrsy = 2e-6; % [m/s]
        end
        tgf0 = min([abs(xG10/vrsy),abs(xG20/vrsy)]);
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
options = odeset("JPattern",jpatt, "RelTol", 1e-2, "AbsTol", 1e-6,...
    "Events",@(t,y)unbond(t,y,ts1,ts2,tp,tc,sTM,sEH,d,Lt,...
                          alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,...
                          b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,...
                          X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,...
                          Fx,Fz,tr,timp,rgf,muT,muG));
%[t,y] = ode45(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0);%, options);
ie = [];
t = [tspan(1), tspan(1)];
y = transpose(y0);
t_out = [];
y_out = [];
t_events = [tspan(1)];
if p.Results.estimateGF
    [t,y,te,ye,ie] = ode45(@(t,y)TMsys_HG(t,y,...
        ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,...
        ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,...
        X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG),...
        [tspan(1),tspan(1)+0.054], y0, options);
    t_out = cat(1, t_out(1:end-1,:), t);
    y_out = cat(1, y_out(1:end-1,:), y);
    xTM = y(:,7:9);
    [yest, ydest] = FitGF(t,xTM,0.004,0.054);
    dfill = (0.01/length(t));
    tfill = ((tspan(1)+0.01):dfill:tgf0)';
    t_out = cat(1, t_out(1:end-1,:), tfill);
    yfill = repmat(y(end,:),length(tfill),1);
    yfill(:,8) = yest(tfill);
    yfill(:,11) = ydest(tfill);
    yfill(:,1) = y(end,1) + vrsy*(tfill-t(end));
    yfill(:,2) = y(end,2) - vrsy*(tfill-t(end));
    y_out = cat(1, y_out(1:end-1,:), yfill);
    y0 = yfill(end,:);
    tspan(1) = tgf0;
    t = [tspan(1), tspan(1)];
else
    y0 = y;
end

while t(end) < tspan(end)
    if numel(ie) > 0
        % Event occurred
        tspan(1) = min(te);
        t_events = cat(1, t_events, tspan(1));
        ie = ie(te<=tspan(1),:);
        ye = ye(te<=tspan(1),:);
        te = te(te<=tspan(1),:);
        y0 = ye(1,:);
        [dydt, dLg1, Fg1, dLg2, Fg2, dLt1, Ft1, dLt2, Ft2] = TMsys_HG(tspan(1),y0',...
            ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,...
            vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,...
            I,T,Fx,Fz,tr,timp,rgf,muT,muG);
        xG1=y_out(end,1);
        xG2=y_out(end,2);
        xT1=y_out(end,3);
        xT2=y_out(end,5);
        xTM=y_out(end,7:9);
        vTM=y_out(end,10:12);
        bTM=y_out(end,13:15);

        fprintf("%f: Release tip contact state is %i, %i.\n", t_out(end), dLt1 < 0, dLt2 < 0);
        
        % No GF contact and 1+ RT separated
        if p.Results.useEnergy && all(ie > 2) && (dLt1 > 0 || dLt2 > 0)
            dT1 = xG1-xT1;
            dT2 = xT2-xG2;
            E = 0.5*kT*dLt1^2*(dLt1<0) + 0.5*kT*dLt2^2*(dLt2<0) + ...
                0.5*k*(dT1^2-d^2)*(dT1>d) + 0.5*k*(dT2^2-d^2)*(dT2>d);
            dv = (-(dLt1<0)+(dLt2<0))*sqrt(2*E/M);
            fprintf("%f: RT impulse = %g.\n", t_out(end), dv);
            maxx = check_capture(xTM, bTM, vTM+[0,dv,0], M, sTM, kGRS);
            captured = sEH/2 - abs(maxx);
            if all(captured>0)
                fprintf("%f: Release tips free and vel within actuation.\n", t_out(end));
                condition = 0;
                break
            end
            if abs(vTM(:,2)) < vrsy
                fprintf("%f: Release tips free and TM vel slower than GF.\n", t_out(end));
                condition = 0;
                break
            end
        end
    else
        if t(1) > tspan(1) && t(end) < tspan(end)
            tspan(1) = t(end);
            y0 = y(1,:);
            fprintf("t = %f: Integration stopped\n",t(end));
            break
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
        if dLt1 < 0
            y0(:,3) = y0(:,3) + dLt1;
            % y0(:,3) = y0(:,8) + sTM/2;
            y0(:,4) = y0(:,11);
            fprintf("%f: RT1/TM in contact\n",te(e5));
        else
            fprintf("%f: RT1/TM lost contact\n",te(e5));
        end
    end
    % RT to TM contact – Set RT vel to TM y
    if any(ie == 6)
        e6 = find(ie==6,1);
        if dLt2 < 0
            y0(:,5) = y0(:,5) - dLt2;
            % y0(:,5) = y0(:,8) - sTM/2;
            y0(:,6) = y0(:,11);
            fprintf("%f: RT2/TM in contact\n",te(e6));
        else
            fprintf("%f: RT2/TM lost contact\n",te(e6));
        end
    end
    % RT to GF contact – Set RT vel to GF, and pos to full length
    if any(ie == 7)
        e7 = find(ie==7,1);
        y0(:,3) = ye(e7,1) - Lt;
        y0(:,4) = vrsy;
        fprintf("%f: GF1/RT1 contact\n",te(e7));
    end
    if any(ie == 8)
        e8 = find(ie==8,1);
        y0(:,5) = Lt - ye(e8,2);
        y0(:,6) = -vrsy;
        fprintf("%f: GF2/RT2 contact\n",te(e8));
    end
    y0(:,23) = y0(:,21) & y0(:,22);
    y0(:,26) = max(y0(:,24),y0(:,25)).*y0(:,23) + -timp.*~y0(:,23);
    [t,y,te,ye,ie] = ode45(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    % [t,y,te,ye,ie] = ode23s(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    % [t,y,te,ye,ie] = ode15s(@(t,y)TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG), tspan, y0, options);
    if numel(ie) > 0 && all(ie < 9) % Don't trim collision events
        t_out = cat(1, t_out(1:end-1,:), t(t<=te(1),:));
        y_out = cat(1, y_out(1:end-1,:), y(t<=te(1),:));
    else
        % No events or contacted EH, so use full series
        t_out = cat(1, t_out(1:end-1,:), t);
        y_out = cat(1, y_out(1:end-1,:), y);
        if numel(ie) > 0
            condition = ie;
            fprintf("Final t=%f.\nEvent t=%f.\n", t(end), te(end));
        else
            condition = 0;
        end
        break
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

maxx = check_capture(xTM, bTM, vTM, M, sTM, kGRS);

%% Output messages

fprintf('\n\nOutput\n\n Max. TM Residual Velocity = %.3d m/s\n',max(abs(vTM(:,2))))
fprintf(' Final TM Residual Velocity = %.3d m/s\n',abs(vTM(end,2)))
if abs(vTM(end,2))>=4.5e-6     % Max residual velocity from "ResVel_2DOF.m"
    fprintf(' Excess y-velocity to avoid collision with EH\n')
elseif condition==9 || abs(maxx(1)) > sEH/2
    fprintf(' TM collision with EH x-face\n')
elseif condition==10 || abs(maxx(2)) > sEH/2
    fprintf(' TM collision with EH y-face\n')
elseif condition==11 || abs(maxx(3)) > sEH/2
    fprintf(' TM collision with EH z-face\n')
else
    fprintf(' TM stopped by actuation with SF = %.1f\n',4.5e-6/abs(vTM(end,2)))
end

% Saved report

% All positive for success
captured = sEH/2 - abs(maxx);

if savedata==1 || saveplot==1 || savestats==1
    mkdir(savepath)
    copyfile(p.Results.config,savepath);
end

if savedata==1
    cd(savepath)
    save(['TMrel_',savingid,'.mat'],'t','xG1','xG2','xT1','vT1','xT2','vT2','xTM','vTM','bTM','wTM','bond1','bond2','-v7.3')
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
    cd('..');
end

if savestats==1
    cd(savepath)
    stats = struct();
    stats.xG1 = get_stats(xG1);
    stats.xG2 = get_stats(xG2);
    stats.xT1 = get_stats(xT1);
    stats.vT1 = get_stats(vT1);
    stats.xT2 = get_stats(xT2);
    stats.vT2 = get_stats(vT2);
    stats.xTM = get_stats(xTM);
    stats.vTM = get_stats(vTM);
    stats.bTM = get_stats(bTM);
    stats.wTM = get_stats(wTM);
    save(['TMrel_',savingid,'_stats.mat'],'stats');
    cd('..');
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
end

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
    cd('..')
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

function [value,isterminal,direction] = unbond(t,y,ts1,ts2,tp,tc,sTM,sEH,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG)
    [dydt, dLg1, Fg1, dLg2, Fg2, dLt1, Ft1, dLt2, Ft2] = TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG);
    bond1 = y(19,1,:);
    bond2 = y(20,1,:);
    % Bonds are broken when exceeding the elongation threshold. To re-form a bond, contact must occur (dLg<=0).
    bond1 = 1*(dLg1<=0)+1*(dLg1>0 & dLg1<tr & bond1==1);
    bond2 = 1*(dLg2<=0)+1*(dLg2>0 & dLg2<tr & bond2==1);
    c1 = (y(1,1,:)-y(3,1,:)-Lt); % GF RT contact
    c2 = (y(5,1,:)-y(2,1,:)-Lt);
    bounds = bbox(y(7:9,1,:), y(13:15,1,:), sTM);
    hs = sEH/2;
    gapx = min([hs-max(bounds(:,1)),hs+min(bounds(:,1))]);
    gapy = min([hs-max(bounds(:,2)),hs+min(bounds(:,2))]);
    gapz = min([hs-max(bounds(:,3)),hs+min(bounds(:,3))]);
    value = [dLg1-tr; dLg2-tr;...   % 1-2
             bond1; bond2;...       % 3-4
             dLt1; dLt2;...         % 5-6
             c1; c2;...             % 7-8
             gapx; gapy; gapz];     % 9-11
    isterminal = [1; 1;...
                  1; 1;...
                  1; 1;...
                  1; 1;...
                  1; 1; 1];
    direction = [1; 1;...
                 0; 0;...
                 0; 0;...
                 1; 1;...
                 -1; -1; -1];
end

function [bounds] = bbox(rTM, bTM, sTM)
    hs = sTM/2;
    corners = [rTM+[hs;hs;hs], rTM+[-hs;hs;hs],...
               rTM+[hs;-hs;hs], rTM+[-hs;-hs;hs],...
               rTM+[hs;hs;-hs], rTM+[-hs;hs;-hs],...
               rTM+[hs;-hs;-hs], rTM+[-hs;-hs;-hs]];
    psi=bTM(1); % angle about x
    the=bTM(2); % angle about y
    phi=bTM(3); % angle about z
    
    % rotational matrix of TM
    R=[cos(the)*cos(phi), cos(the)*sin(phi) ,-sin(the);...
       sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi),...
       sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi),...
       sin(psi)*cos(the);...
       cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi),...
       cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi),...
       cos(psi)*cos(the)];
    bounds = (R*corners)';
end

function [maxx] = check_capture(xTM, bTM, vTM, M, sTM, kGRS)
    bounds = bbox(xTM', bTM', sTM);
    edges = [min(bounds(:,1)), min(bounds(:,2)), min(bounds(:,3));...
             max(bounds(:,1)), max(bounds(:,2)), max(bounds(:,3))];
    posv = (vTM(end,:)>0);
    % Find extent of sides in direction of travel
    lead_sides = edges(1,:);
    lead_sides(posv) = edges(2,posv);
    maxarg = vTM(end,:)./lead_sides .* sqrt(M/abs(kGRS));
    maxx = lead_sides./sqrt(1+maxarg.^2) + vTM(end,:).*sqrt(M/abs(kGRS)).*maxarg./sqrt(1+maxarg.^2);
end

function [stats] = get_stats(var)
    stats = struct();
    stats.max = max(abs(var));
    stats.final = var(end,:);
    stats.mean = mean(var);
    stats.std = std(var);
end