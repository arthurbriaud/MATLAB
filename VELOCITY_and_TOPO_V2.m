close all;
clear;
clc;

%--------------------------------------------------------------------------------------------%
% Biggest code divided in many parts where topography,isostasy,dynamic topography, ...
% locations and velocities are computed.
/%--------------------------------------------------------------------------------------------%

main_dir = '/Users/Arthur workstation/Documents/SIMULATIONS/';

%---------------------------
% Step 
%---------------------------
Step_init = 0;
Step_Step = 100;
Step_end  = 0;
%---------------------------

write = 0;
fig = 0;

pro = 96;k = 1e-6; DT = 1300; H = 3000e3; g = 9.8;ro0=3300; alfa = 3e-5;
mu_ref = 1e20; dx=1; roM= 3300; row = 1000; roC = 2700; 
Ridge_alt = 2.9;  %[km]  


fid = fopen('C:\Users\Arthur workstation\Documents\MATLABpost_pro\Data\NaN.txt','w');

for i =  Step_init : Step_Step : Step_end
disp('-------------------------');
disp([' Plotting Step nb: ',num2str(i,'%04d')]);
disp('-v-v-v-v-v-v-v-v-v-v-v-v-');
 cd(strcat(main_dir, 'MODEL NAME'))
 
    [V,P,Surf,C,t] = CITCOM_reader(i,pro,1,2);
    X = C{1} .*H.*1e-3; 
    Z = C{2} .*H.*1e-3;
    T = V{1}.*DT; 
    mu = V{4}.*mu_ref;
    t = t{1}.*(H^2/k) ./ 31536000e6;   % [Myr]
    Vx = V{2}.*(k/H) .* 31536000e2; % [cm/yr]
    Vz = V{3}.*(k/H) .* 31536000e2; % [cm/yr]
    Vxx = Vx .*-1; % Vx positif towards left
    Rheol = V{5};
    Vx2 = Vx.^2; Vz2 = Vz.^2;  % squared   % [cm/yr]
    V_rms = sqrt(Vx2 + Vz2);    % Vrms magnitude % [cm/yr]
    
   %%  | Topography/Stress |     
    Topo_oc = Surf{1}.*(mu_ref*k/H^2)./(roM-row)/g.*1e-3;     % [km] considering water above the oceaninc lithosphere
    Topo_cc = Surf{1}.*(mu_ref*k/H^2)./(roC-1)/g.*1e-3;       % [km] maybe consider air above the continent!

    % find trench with rheology is easier (just look for V{5}==0 (above subd ptale)
    % or V{5}==-2 above upper plate) 
    tr_idx = find(Rheol(1,:)==0,1,'last');
    % tr_idx_T = find(Rheol(1,:)==-2,1,'first');
    TOPO_oc = Topo_oc(1:tr_idx,1);          % Upper plate topography [km]
    TOPO_cc = Topo_cc(tr_idx+1:end,1);      % Lower plate topography [km]

    TOPO = [TOPO_oc;TOPO_cc];   % Actual topography
    
 % scale topo regarding to the actual ridge depth
    TOPO_cor = TOPO(1:tr_idx-200,1);
    X_oc = X(1:tr_idx-200,1);
    pmax = max(TOPO_cor);
    [~,id] = findpeaks(TOPO_cor,'MINPEAKHEIGHT',0.8*pmax);
    cor1 = pmax + Ridge_alt;
    Topo_OP = TOPO - cor1; 
    TOPO    = TOPO - cor1;
    
    %% | Isostatic balance and Dynamic topography |

    idx_comp = find(Z>=300,1,'first');  % find index for compesation depth of 300 km
    TOC = T(1:idx_comp,1:tr_idx);       % temperaure field for the oceanc side
    TC = T(1:idx_comp,tr_idx+1:end);    % temperature filed for the continental side
    TM = T(1:idx_comp,tr_idx+1:end);     % mantle temperature field
    RH = Rheol(1:idx_comp,tr_idx+1:end);
    TC(RH~=-2) = 0; TC(RH==-2) = 1;

    % for the oceanic side, the stress is:
%     sigmaZZ = int(rho_mantle * g * (alfa ( T oceanic - DT))
    Soc = trapz(Z(1:idx_comp).*1e3,g.*roM.*(alfa.*(TOC-DT)));
    % for the continent side the stress is:
%     sigmaZZ = rho_mantle*g*H_all + (roM-roC)*g*H_continent
    Scc = trapz(Z(1:idx_comp).*1e3,g.*(roM-roC).*TC) + ...
          trapz(Z(1:idx_comp).*1e3,g.*roM.*(alfa.*(TM-DT)));
    
    % dynamic topography is surface topography - isostasy 
    TOPoc_dyn = TOPO_oc - Soc'./(roM-row)/g.*1e-3;
    TOPcc_dyn = TOPO_cc - Scc'./(roC-1)/g.*1e-3;
    TOPO_dyn  = [TOPoc_dyn;TOPcc_dyn];

%% Find locations 
% Trench location
Mu_tmp = mu(1:1,1:end);      
[~,tr_idx_tmp]     = find(Mu_tmp>=2e20);        % 5e19 [Pa*s] weak zone reference for oceanic upper plate % 2e20 for continental trench location
tr_idx_1             = tr_idx_tmp(1);           % [node]  (1) if continental crust is here 
Tre_x              = X(tr_idx_1);               % [km]
Edge_x_tmp         = tr_idx_tmp(end);           % [node]
Edge_x             = X(Edge_x_tmp);             % [km]
tr_idx_T1          = 2248;                      % first trench location
    
x                   = [1:dx:X(end)];
Mu_tmp2             = interp1(X,Mu_tmp,x);
[~,tr_idx_tmp2]     = find(Mu_tmp2>=2e20);
tr_idx2             = tr_idx_tmp2(1);                  % [km]
Topo_OP_tmp         = interp1(X,Topo_OP,x);            % [km]
Topo_OP             = Topo_OP_tmp(1,tr_idx2:end);      % [km]
Xtrench_topo        = Topo_OP_tmp(1,tr_idx2);          % [km]
TOPO_dyn_tmp        = interp1(X,TOPO_dyn,x);
% Upper plate location rear the trench
OVPref = [tr_idx2+800 tr_idx2+2200] ; 
OVP =    x(OVPref(1):OVPref(2));

% Subducting plate location before the trench
SVP_ref = [tr_idx2-800];
%% Averaged velocities and surface processes 

V_tmp        = Vx(1:1,1:end);
Vx_tmp2      = interp1(X,V_tmp,x); 
V_OP         = mean(Vx_tmp2(1,OVPref(1):OVPref(2))); % Upper plate velocity [cm/yr]
V_trench     = mean(Vx_tmp2(tr_idx2));           % Trench velocity [cm/yr]
V_sub        = Vx_tmp2(SVP_ref);                 % subducting plate velocity [cm/yr]
dh_surf_topo = (Topo_OP_tmp(OVPref(2)) - Topo_OP_tmp(OVPref(1)));      % [km]

dX          = (x(OVPref(2)) - x(OVPref(1)));                           % [km]
Tilt_OP     = dh_surf_topo/dX;               % surface topograhy tilting, slope [-]
Altitude_OP = mean(Topo_OP_tmp(OVP));        % surface topography up and down [km] ;
TOPO_dyn_OP = mean(TOPO_dyn_tmp(OVP));       % averaged dynamic continental topography up and down [km]; 
disp('TOPO_dyn_OP')
disp(TOPO_dyn_OP)

 %% Define lithosphere location 
T_refVrms = zeros(size(T));
  for m=1:numel(Z)
           for n=1:numel(X)
                 T_ref_tmp = (T(m,n)>=1200); % temperature higher than 1200 C  
                 T_refVrms(m,n)= T_ref_tmp; 
            end
  end
clear T_ref_tmp   
%% Define only the slab
T_refVs = zeros(size(T));
  for m=1:numel(Z)
           for n=1:numel(X)
                 T_ref_tmp = (T(m,n)<=1200); % temperature lower than 1200 C  
                 T_refVs(m,n)= T_ref_tmp; 
            end
  end
clear T_ref_tmp
T_refVr = T_refVs;
%% 
idx = find(Z>=0,1,'first'); % vertical index where Z>= 0 km
idx1 = find(Z>=200,1,'first'); % vertical index where Z>= 200 km
idx2 = find(Z>=670,1,'first'); % vertical index where Z>= 670 km

%% Vrms

Z_upper = Z(idx:idx2);
Z_lower = Z(idx2:1:end);

V_RMS = V_rms;        
V_RMS(T_refVrms==0) = NaN; 
T_refVrms(T_refVrms==0) = NaN;  
V_RMS_upper = V_RMS(idx:idx2,:); % Vrms in the upper mantle
V_RMS_lower = V_RMS(idx2:end,:);  % Vrms in the lower mantle

V_RMS_upper(isnan(V_RMS_upper)) = 0;
V_RMS_lower(isnan(V_RMS_lower)) = 0;

V_RMS_upper = V_RMS_upper * 3.1710e-6;   % [m/s]
V_RMS_lower = V_RMS_lower * 3.1710e-6;   % [m/s]

INT_upper_tmp = trapz(Z_upper.*1e3,V_RMS_upper);
INT_lower_tmp = trapz(Z_lower.*1e3,V_RMS_lower);
INT_upper = INT_upper_tmp.'; INT_lower = INT_lower_tmp.'; 

%% Sinking velocity
VZ = Vz;   
VZ(T_refVs==0) = NaN; 
T_refVs(T_refVs==0) = NaN;  
VZ_upper = VZ(idx1:idx2,:);
VZ_lower = VZ(idx2:end,:); Z_roll = Z(idx1:idx2,:);

%% Roll Back velocity
VX = Vx;
VX(T_refVr==0) = NaN;
T_refVr(T_refVr==0) = NaN;
V_ROLL_back = VX(idx1:idx2,:);

%% calculate averaged Vrms velocities
Vrms_upper = V_RMS_upper(isfinite(V_RMS_upper)); % retourn only the finite value
Vrms_lower = V_RMS_lower(isfinite(V_RMS_lower));
Vrms_Upper = mean(Vrms_upper);   
Vrms_Lower = mean(Vrms_lower);   

 % coupling decoupling between both layered mantle 
 
Vrms_up = mean(mean(INT_upper));    % integrate velocity on the upper mantle area
Vrms_dm = mean(INT_lower);          % for the lower mantle
Coupling = Vrms_up/Vrms_dm;         % Coupling/decoupling between upper and lower mantle 

disp('Vrms Upper')
disp(Vrms_Upper)
disp('Vrms_Lower')
disp(Vrms_Lower)
disp('Coupling Decoupling')
disp(Coupling)
figure(1)
subplot(3,3,[1 3])
plot(t,Coupling,'rs'); title('Coupling'); hold on; 
subplot(3,3,[4 6])
plot(t,Altitude_OP,'kv');title('Alt'); hold on;
subplot(3,3,[7 9])
plot(t,Tilt_OP,'bv');title('Tilt'); hold on;
drawnow

%% calculate averaged sinking velocities

Vsink_upper = VZ_upper(isfinite(VZ_upper)); % retourn only the finite value
Vsink_lower = VZ_lower(isfinite(VZ_lower));

Vsink_Upper = mean(Vsink_upper);    % note: the value will be NaN if the matric is empity, 
Vsink_Lower = mean(Vsink_lower);    %       in other words if the slab is not present at these depths

disp('Vsink_Upper')
disp(Vsink_Upper)
disp('Vsink_Lower')
disp(Vsink_Lower)

figure(2)
subplot(3,3,[1 3])
plot(t,Vsink_Upper,'rs',t,Vsink_Lower,'bs'); title('Sinking'); hold on; 
subplot(3,3,[4 6])
plot(t,TOPO_dyn_OP,'kv');title('dyn top alt'); hold on;
subplot(3,3,[7 9])
plot(t,Tilt_OP,'bv');title('Tilt'); hold on;
drawnow

%% calculate averaged roll back velocity

V_roll_back = V_ROLL_back(isfinite(V_ROLL_back)); % retourn only the finite value
V_roll_back = mean(V_roll_back); 
disp('Roll back')
disp(V_roll_back)

figure(3)
subplot(2,2, [1 2])
plot(t,V_roll_back,'ms',t,V_OP,'ks',t,V_sub,'bs'); title('Roll back (m) & upper plate (k) & V sub'); hold on; 
subplot(2,2,[3 4]);
plot(t,-V_OP,'ks'); title('V OP'); hold on;
plot(t,-V_trench,'bv'); title('V trench'); hold on;
%% Figure
if (fig ==1)
%% Vrms

%% first check : Vrms only within the mantle 

figure(i+1)
subplot(221)
pcolor(X,Z,V_rms); shading interp; colormap jet; axis ij; axis equal;
colormap jet; colorbar; hold on; %caxis([0 8]);
title('Vrms all system')
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(222),pcolor(X,Z,V_RMS); axis ij; shading interp; axis equal;
colormap jet; colorbar; hold on; %caxis([0 8]);
title('V rms')
for(pp=1:2),contour(X,Z,P{pp},[-0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(223),pcolor(X,Z(idx:idx2),V_RMS_upper); axis ij; shading interp; axis equal;
colormap jet; colorbar; hold on; %caxis([0 8]);
title('V rms upper')
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(224),pcolor(X,Z(idx2:end),V_RMS_lower); axis ij; shading interp;axis equal;
colormap jet; colorbar; hold on; % caxis([0 8]);
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end
title('V rms lower')
drawnow;

%% Sinking speed 

%% Sinking speed

figure(i+10000) %  Isotherm profile
subplot(221),pcolor(X,Z,T_refVs); axis ij; shading interp;
colormap jet; colorbar; axis equal; hold on
title('T ref');
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(222),pcolor(X,Z,VZ); axis ij; shading interp;
colormap jet; colorbar; axis equal; hold on;
title('VZ')
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(223),pcolor(X,Z(idx1:idx2),VZ_upper); axis ij; shading interp;
colormap jet; colorbar; axis equal; hold on;
title('VZ upper')
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end

subplot(224),pcolor(X,Z(idx2:end),VZ_lower); axis ij; shading interp;
colormap jet; colorbar; axis equal; hold on
for(pp=1:2),contour(X,Z,P{pp},[0.5,0.5],'k','LineWidth',1); axis ij; end
title('VZ lower')
drawnow;
end


      if(write ==1)
 fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',t,Coupling,Vrms_up,Vrms_dm,Vsink_Upper,Vsink_Lower,Altitude_OP,Tilt_OP,Tre_x,V_roll_back,V_OP,V_trench,V_sub,i,TOPO_dyn_OP); 
      end
      
      
end
