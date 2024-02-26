 %% PSSL - ILaRIS CVM
% Function that describes Actuator, Grabbing Finger (GF), Release Tip (RT), and Test Mass (TM) system
% according to paper "Prediction of the LISA-Pathfinder release mechanism
% in-flight performance" and the spring-tip CVM concept. Called by the script "TMrel_Ball3.m"

function [dydt, dLg1, Fg1, dLg2, Fg2, dLt1, Ft1, dLt2, Ft2]= TMsys_HG(t,y,ts1,ts2,tp,tc,sTM,d,Lt,alp,ac,bc,mT,M,k,kG,kT,kGRS,ixp,ixn,izp,izn,iy,b,vrsx,vrsy,X1g1,X2g1,X3g1,X1g2,X2g2,X3g2,X1t1,X2t1,X3t1,X1t2,X2t2,X3t2,I,T,Fx,Fz,tr,timp,rgf,muT,muG,estimateRT)
% global imp imn imy bond1 bond2 tsep1 tsep2 tsept

nvar = numel(y);
bond1 = y(19);
bond2 = y(20);
imp = y(21);
imn = y(22);
imy = y(23);
tsep1 = y(24);
tsep2 = y(25);
tsept = y(26);

%% Input
% GF speed - assuming constant actuator speed

% Find value in time-series
t1 = t-ts1; % time since delay (to use GF speed profiles at proper times)
t2 = t-ts2;

% Option 1 - interpolate the current speed from an experimentally obtained time-series
% vr1f = interp1(vrfx,vrfy,t1,'nearest','extrap'); % fast and slow GF1 speed at time
% % vr1s = interp1(vrsx,vrsy,t1,'nearest','extrap');
% vr2f = -1*interp1(vrfx,vrfy,t2,'nearest','extrap'); % fast and slow GF2 speed at time
% vr2s = interp1(vrsx,vrsy,t2,'nearest','extrap');
% vG1 = (t>=ts1)*vr1s; % current GF speeds
% vG2 = -1*(t>=ts2)*vr2s;

% Option 2 - start with a speed and change to another
% vG1 = vrfy*(t>=ts1 & t<tp)+vrsy*(t>=tp); % current GF speeds
% vG2 = -1*(vrfy*(t>=ts2 & t<tp)+vrsy*(t>=tp));

% Option 3 - use a constant speed value
vG1 = (t>=ts1)*vrsy; % current GF speeds
vG2 = -1*(t>=ts2)*vrsy;

%% Adjusting offset (comment out if 15 mm offset is added as initial condition)
yo = zeros([5,1]);
yo(1)=y(1)+sTM/2;
yo(2)=y(2)-sTM/2;
yo(3)=y(3)+sTM/2;
yo(5)=y(5)-sTM/2;

% yo(1)=y(1);
% yo(2)=y(2);
% yo(3)=y(3);
% yo(5)=y(5);

%% TM Orientation

psi=y(13); % angle about x
the=y(14); % angle about y
phi=y(15); % angle about z

% rotational matrix of TM
R=[cos(the)*cos(phi) , cos(the)*sin(phi) , -sin(the) ; sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi) , sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) , sin(psi)*cos(the) ; cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi) , cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) , cos(psi)*cos(the)];

s1 = (R*[0,sTM/2,0]')'; % vector to +Y face
s2 = (R*[0,-sTM/2,0]')'; % vector to -Y face
s1h = s1/norm(s1); % unit vector of s1
s2h = s2/norm(s2); % unit vector of s2

%% Position and Speed Vectors

yh = [0 1 0]; % Y-axis

rt1 = [0 yo(3) 0]; % RT1 position
rt2 = [0 yo(5) 0]; % RT2 position
vt1 = [0 y(4) 0]; % RT1 speed
vt2 = [0 y(6) 0]; % RT2 speed

rg1 = [0 yo(1) 0]; % GF1 position
rg2 = [0 yo(2) 0]; % GF2 position
vg1 = [0 vG1 0]; % GF1 speed
vg2 = [0 vG2 0]; % GF2 speed

rTM = [y(7),y(8),y(9)]; % TM position
omTM = [y(16),y(17),y(18)]; % TM angular speed [rad/s]
vTM = [y(10),y(11),y(12)]; % TM linear speed [m/s]

%% Contact Forces

% GF to RT contact
Fc1=kG*(yo(1)-Lt-yo(3))*(yo(1)-yo(3)>Lt);
Fc2=kG*(yo(5)-Lt-yo(2))*(yo(5)-yo(2)>Lt);

% GF to TM contact
ys1 = (yh-(yh*s1h')*s1h)/norm(yh-(yh*s1h')*s1h); % Y-axis component on TM Y-face
ys2 = -1*((yh-(yh*s2h')*s2h)/norm(yh-(yh*s2h')*s2h));

% Contact point vector (comment to assume contact occurs at GF center)
% rp1 = rgf*((ys1-(ys1*yh')*yh)/norm(ys1-(ys1*yh')*yh)); % Contact point from GF
% rp2 = rgf*((ys2-(ys2*yh')*yh)/norm(ys2-(ys2*yh')*yh));
% rp1(isnan(rp1))=0; % Checking for NaNs (occurs if Y-face normal is parellel to fixed Y-axis)
% rp2(isnan(rp2))=0;
% lg1 = rp1+rg1-(rTM+s1); % GF contact point from TM Y-face center
% lg2 = rp2+rg2-(rTM+s2);

lg1 = rg1-(rTM+s1); % GF contact point from TM Y-face center
lg2 = rg2-(rTM+s2);
dLg1=lg1*s1h'; % GF to TM face distance along TM face normal
dLg2=lg2*s2h';

% Bonds are broken when exceeding the elongation threshold. To re-form a bond, contact must occur (dLg<=0).
bond1 = 1*(dLg1<=0)+1*(dLg1>0 && dLg1<tr && bond1==1);
bond2 = 1*(dLg2<=0)+1*(dLg2>0 && dLg2<tr && bond2==1);

Fg1 = kG*dLg1*(dLg1<=0) + X1g1*dLg1.*exp(-X2g1*dLg1.^X3g1)*(dLg1>0 && bond1==1);
Fg2 = kG*dLg2*(dLg2<=0) + X1g2*dLg2.*exp(-X2g2*dLg2.^X3g2)*(dLg2>0 && bond2==1);
% Fg1 = (kG*dLg1*(dLg1<=0) + X1g1*dLg1.*exp(-X2g1*dLg1.^X3g1)*(dLg1>0 && bond1==1))*(yo(5)-yo(2)>Lt-5e-6); % MOD
% Fg2 = (kG*dLg2*(dLg2<=0) + X1g2*dLg2.*exp(-X2g2*dLg2.^X3g2)*(dLg2>0 && bond2==1))*(yo(5)-yo(2)>Lt-5e-6); % MOD

% RT to TM contact
lt1 = rt1-rTM-s1; % RT contact point from TM Y-face center
lt2 = rt2-rTM-s2;
dLt1 = lt1*s1h'; % RT to TM face distance along TM face normal
dLt2 = lt2*s2h';

% Check case conditions % MOD
cT1 = dLt1<0; % RT1 to TM contact
cT2 = dLt2<0; % RT2 to TM contact
cG1 = Fc1>0; % GF1 to RT1 contact
cG2 = Fc2>0; % GF2 to RT2 contact

Ft1=kT*dLt1*(dLt1<=0) + X1t1*dLt1.*exp(-X2t1*dLt1.^X3t1)*(dLt1>0); 
Ft2=kT*dLt2*(dLt2<=0) + X1t2*dLt2.*exp(-X2t2*dLt2.^X3t2)*(dLt2>0);
% Ft1=(kT*dLt1*(dLt1<=0) + X1t1*dLt1.*exp(-X2t1*dLt1.^X3t1)*(dLt1>0))*(cT1<=cG1) + k*d*(cT1>cG1); % MOD default eq. except when there is only contact between RT and TM
% Ft2=(kT*dLt2*(dLt2<=0) + X1t2*dLt2.*exp(-X2t2*dLt2.^X3t2)*(dLt2>0))*(cT2<=cG2) + k*d*(cT2>cG2); % MOD


%% Friction Direction Vectors

lts1 = lt1-(lt1*s1h')*s1h; % vector along Y-face to contact point
lts2 = lt2-(lt2*s2h')*s2h;
lgs1 = lg1-(lg1*s1h')*s1h;
lgs2 = lg2-(lg2*s2h')*s2h;

rtc1 = s1+lts1; % vector from TM center to contact point
rtc2 = s2+lts2;
rgc1 = s1+lgs1;
rgc2 = s2+lgs2;

vtc1 = cross(omTM,rtc1)+vTM; % absolute speed of contact point
vtc2 = cross(omTM,rtc2)+vTM;
vgc1 = cross(omTM,rgc1)+vTM;
vgc2 = cross(omTM,rgc2)+vTM;

mt1v = (vtc1-vt1)/norm(vtc1-vt1); % relative speed unit vector of contact point to colliding object (GF or RT)
mt2v = (vtc2-vt2)/norm(vtc2-vt2);
mg1v = (vgc1-vg1)/norm(vgc1-vg1);
mg2v = (vgc2-vg2)/norm(vgc2-vg2);

mt1 = (mt1v-(mt1v*s1h')*s1h)/norm(mt1v-(mt1v*s1h')*s1h); % relative speed unit vector tangent to TM Y faces
mt2 = (mt2v-(mt2v*s2h')*s2h)/norm(mt2v-(mt2v*s2h')*s2h);
mg1 = (mg1v-(mg1v*s1h')*s1h)/norm(mg1v-(mg1v*s1h')*s1h);
mg2 = (mg2v-(mg2v*s2h')*s2h)/norm(mg2v-(mg2v*s2h')*s2h);

mt1(isnan(mt1))=0; % Checking for NaNs (occurs if the relative speed of the TM to the colliding object is 0)
mt2(isnan(mt2))=0;
mg1(isnan(mg1))=0;
mg2(isnan(mg2))=0;

%% TM Force Vectors

% Normal forces [N]
Fg1v = Fg1*s1h;
Fg2v = Fg2*s2h;
Ft1v = Ft1*s1h;
Ft2v = Ft2*s2h;

% Tangential forces [N]
Ffg1 = -muG*abs(Fg1)*mg1;
Ffg2 = -muG*abs(Fg2)*mg2;
Fft1 = -muT*abs(Ft1)*mt1;
Fft2 = -muT*abs(Ft2)*mt2;

Ftm = Fg1v+Ft1v+Fg2v+Ft2v+Ffg1+Fft1+Ffg2+Fft2; % Sum of all contact forces on TM
Ttm = -1*(cross(lgs1+s1,Fg1v)-cross(lgs1+s1,Ffg1)+cross(lts1+s1,Ft1v)-cross(lts1+s1,Fft1)+cross(lgs2+s2,Fg2v)-cross(lgs2+s2,Ffg2)+cross(lts2+s2,Ft2v)-cross(lts2+s2,Fft2)); % Sum of all contact torques on TM

%% Lateral Impulse Time
% Lateral impulses are applied when the GF loses contact with the TM for a
% single step (based on in-flight LPF GPRM data)

if dLg1>tr && imp==0
    tsep1 = t; % time of separation of GF1
    imp = 1;
end
% tsep1 = t*(dLg1>tr && imp==0) -timp*~(dLg1>tr && imp==0);
% imp = 1*(dLg1>tr && imp==0);

if dLg2>tr && imn==0
    tsep2 = t; % time of separation of GF2
    imn = 1;
end

if imp==1 && imn==1 && imy==0
    tsept = t; % time of separation from both GFs
    imy = 1;
end

% CGF side (GF1) impulse
ixp1 = ixp*(t-tsep1<timp);
izp1 = izp*(t-tsep1<timp);

% PGF side (GF2) impulse
ixn1 = ixn*(t-tsep2<timp);
izn1 = izn*(t-tsep2<timp);

% Y-axis impulse (once both GF1 and GF2 are out of contact)
iy1=iy*(t-tsept<timp);

%% System Equations
% state vector y = [xG1 xG2 xT1 vT1 xT2 vT2 xTM yTM zTM v_xTM v_yTM v_zTM an_xTM an_yTM an_zTM om_xTM om_yTM om_zTM]

% Grabbing Finger and Release Tip (1 DOF, Y-axis displacement)
    dydt = zeros([nvar,1]);
    estm = 10^(-0.99096642*log10(kT)-1.640234);
    % estb = 10^(-0.98271567*log10(kT)-0.41746246);
    % estR1 = t*estm-estb
    % estR2 = -t*estm+estb
    dydt(1,:) = vG1;                                                               % vG1
    dydt(2,:) = vG2;                                                               % vG2
    
    % if cT1
    %     dydt(3,:) = y(11);
    %     if y(4) > 0
    %         dydt(3,:) = 0;
    %     end
    %     dydt(4,:) = 0;
    % else
        dydt(3,:) = y(4); % vT1
    if cG1
        dydt(4,:) = 0;
    elseif estimateRT > 0 && cT1(1) && y(3) < -0.25e-7
        dydt(3,:) = estm;
        dydt(4,:) = 0;
    else
        dydt(4,:) = (1/mT)*(k*(yo(1)-yo(3)-d)+b*(vG1-y(4))+Fc1-Ft1v(2)); % aT1
    end
    % end
    % if cT2
    %     dydt(5,:) = y(11);
    %     if y(6) > 0
    %         dydt(5,:) = 0;
    %     end
    %     dydt(6,:) = 0;
    % else
        dydt(5,:) = y(6); % vT2
    if cG2
        dydt(6,:) = 0;
    elseif estimateRT > 0 && cT2(1) && y(5) > 0.25e-7
        dydt(5,:) = -estm;
        dydt(6,:) = 0;
    else
        dydt(6,:) = (1/mT)*(k*(yo(2)-yo(5)+d)+b*(vG2-y(6))-Fc2-Ft2v(2)); % aT2
    end
    % end

  % Test Mass (6 DOF)
    dydt(7,:) = y(10);                                      % vxTM
    dydt(8,:) = y(11);                                                % vyTM
    dydt(9,:) = y(12);                                        % vzTM
    dydt(10,:) = (1/M)*(Ftm(1)+Fx-kGRS*y(7)+(ixp1+ixn1)/timp); % axTM
    dydt(11,:) = (1/M)*(Ftm(2)-kGRS*y(8)+iy1/timp); % ayTM
    dydt(12,:) = (1/M)*(Ftm(3)+Fz-kGRS*y(9)+(izp1+izn1)/timp); % azTM
    dydt(13,:) = y(16); % wxTM
    dydt(14,:) = y(17); % wyTM
    dydt(15,:) = y(18); % wzTM
    dydt(16,:) = (1/I(1,1))*(Ttm(1)+T(1)+(-izp1+izn1)*(ac+bc*tan(alp))/timp)+(I(2,2)-I(3,3))*y(17)*y(18); % alxTM
    dydt(17,:) = (1/I(2,2))*(Ttm(2)+T(2))+(I(3,3)-I(1,1))*y(16)*y(18);                                      % alyTM
    dydt(18,:) = (1/I(3,3))*(Ttm(3)+T(3)+(ixp1-ixn1)*(ac+bc*tan(alp))/timp)+(I(1,1)-I(2,2))*y(16)*y(17);  % alzTM

    % if t > 2.6% && mod(floor(t*100),10)==0
    %     plot(rtc1(1), rtc1(3), ".b");
    %     hold on;
    %     % plot(t, mt1v(2), ".r");
    %     % plot(t, mt1v(3), ".g");
    %     % xlim([-1,t]);
    %     % ylim([-.1,.1]);
    %     drawnow;
    % end
    
% MODIFIED EoMs (speeds code with contact assumptions, inaccurate)
%     dydt(3,:) = y(4)*(cT1==cG1) + y(11)*cT1*~cG1 + vG1*cG1*~cT1;                   % vT1  
%     dydt(4,:) = (1/mT)*(k*(yo(1)-yo(3)-d)+b*(vG1-y(4))+Fc1-Ft1v(2))*(cT1==cG1);    % aT1
%     dydt(5,:) = y(6)*(cT2==cG2) + y(11)*cT2*~cG2 + vG2*cG2*~cT2;                   % vT2
%     dydt(6,:) = (1/mT)*(k*(yo(2)-yo(5)+d)+b*(vG2-y(6))-Fc2-Ft2v(2))*(cT2==cG2);    % aT2
%     dydt(7,:) = y(10)*(yo(5)-yo(2)>Lt-5e-6);                                      % vxTM
%     dydt(8,:) = y(11)*(yo(5)-yo(2)>Lt-5e-6)+vG2*(yo(5)-yo(2)<=Lt-5e-6);                                                % vyTM
%     dydt(9,:) = y(12)*(yo(5)-yo(2)>Lt-5e-6);                                        % vzTM 
%     dydt(10,:) = (1/M)*(Ftm(1)+Fx-kGRS*y(7)+(ixp1+ixn1)/timp)*(yo(5)-yo(2)>Lt-5e-6); % axTM 
%     dydt(11,:) = (1/M)*(Ftm(2)-kGRS*y(8)+iy1/timp)*(yo(5)-yo(2)>Lt-5e-6); % ayTM 
%     dydt(12,:) = (1/M)*(Ftm(3)+Fz-kGRS*y(9)+(izp1+izn1)/timp)*(yo(5)-yo(2)>Lt-5e-6); % azTM 
%     dydt(13,:) = y(16)*(yo(5)-yo(2)>Lt-5e-6); % wxTM
%     dydt(14,:) = y(17)*(yo(5)-yo(2)>Lt-5e-6); % wyTM
%     dydt(15,:) = y(18)*(yo(5)-yo(2)>Lt-5e-6); % wzTM
%     dydt(16,:) = ((1/I(1,1))*(Ttm(1)+T(1)+(-izp1+izn1)*(ac+bc*tan(alp))/timp)+(I(2,2)-I(3,3))*y(17)*y(18))*(yo(5)-yo(2)>Lt-5e-6); % alxTM
%     dydt(17,:) = ((1/I(2,2))*(Ttm(2)+T(2))+(I(3,3)-I(1,1))*y(16)*y(18))*(yo(5)-yo(2)>Lt-5e-6);                                      % alyTM
%     dydt(18,:) = ((1/I(3,3))*(Ttm(3)+T(3)+(ixp1-ixn1)*(ac+bc*tan(alp))/timp)+(I(1,1)-I(2,2))*y(16)*y(17))*(yo(5)-yo(2)>Lt-5e-6);  % alzTM
    
%% DEBUG
% dst = 300; 
% plot(tout(1:end-dst),yout(1,1:end-dst),tout(1:end-dst),yout(2,1:end-dst),tout(1:end-dst),yout(3,1:end-dst),tout(1:end-dst),yout(5,1:end-dst),tout(1:end-dst),yout(8,1:end-dst))
% legend GF1 GF2 RT1 RT2 TM