close all; clear; clc;
%-------------------------------------------------------------------------------------%
% this code is used for trench retreat and advance location according to trench velocities.
%-------------------------------------------------------------------------------------%

% load file
load cC50cO100g6_6.txt
t_1_1             = cC50cO100g6_6(1:end,1); % [Myr]
Coupling_1_1      = cC50cO100g6_6(1:end,2); % [-]
Vrms_Upper_1_1    = cC50cO100g6_6(1:end,3); % [cm/yr]
Vrms_Lower_1_1    = cC50cO100g6_6(1:end,4); % [cm/yr]
Vsink_Upper_1_1   = cC50cO100g6_6(1:end,5); % [cm/yr]
Vsink_Lower_1_1   = cC50cO100g6_6(1:end,6); % [cm/yr]
Altitude_OP_1_1   = cC50cO100g6_6(1:end,7); % [km]
Tilt_OP_1_1       = cC50cO100g6_6(1:end,8); % [km/km]
Tr_location_1_1   = cC50cO100g6_6(1:end,9); % [km]
V_rollback_1_1    = cC50cO100g6_6(1:end,10);% [cm/yr]
V_upperplate_1_1  = cC50cO100g6_6(1:end,11);% [cm/yr]
V_trench_1_1      = cC50cO100g6_6(1:end,12);% [cm/yr]
V_sub_1_1         = cC50cO100g6_6(1:end,13);% [cm/yr]
Ii                 = cC50cO100g6_6(1:end,14); %[iterations]
Topo_dyn_1_1      = cC50cO100g6_6(1:end,15); %[km]
TAB = table(t_1_1,Ii);

% trench location nomalised to the first step wanted
trench1_1 = Tr_location_1_1(1:1)-Tr_location_1_1; % trench position at t=0 Myr; x = 0km
XLoc1_1  = trench1_1(end,end);

% get trench retreat and advance
% retreat
tr_rt_1_1 = trench1_1; t_rt_1_1 = t_1_1; % save time and trench vectors
Vtr_rt_1_1 = -V_upperplate_1_1; % define trench velocity i.e., same tham V_trench
tr_rt_1_1(Vtr_rt_1_1<=0) = NaN; % NaN when trench velocity is negative for trench location
Vtr_rt_1_1(Vtr_rt_1_1<=0) = NaN;% NaN when trench velocity is negative for trench velocity
t_rt_1_1(isnan(tr_rt_1_1)) = NaN; % NaN when trench velocity is negative for time
% advance
tr_ad_1_1 = trench1_1; t_ad_1_1 = t_1_1; % save time and trench vectors
Vtr_ad_1_1 = -V_upperplate_1_1; % define trench velocity i.e., same tham V_trench
tr_ad_1_1(Vtr_ad_1_1>=0) = NaN;% NaN when trench velocity is positive for trench location
Vtr_ad_1_1(Vtr_ad_1_1>=0) = NaN; % NaN when trench velocity is positive for trench velocity
t_ad_1_1(isnan(tr_ad_1_1)) = NaN; % NaN when trench velocity is positive for time

% check trench retreat and advance velocities
    % blue = retreat red = advance
figure(100)
plot(t_rt_1_1,Vtr_rt_1_1,'bs'); hold on; if(numel(Vtr_ad_1_1~=0)); plot(t_ad_1_1,Vtr_ad_1_1,'rs'); end; hold on;
plot(t_1_1,-V_upperplate_1_1,'k');
% check trench displacment too
figure(101)
plot(t_rt_1_1,tr_rt_1_1,'bs'); hold on; if(numel(Vtr_ad_1_1~=0)); plot(t_ad_1_1,tr_ad_1_1,'rs'); end

tr_rt_1_1 = tr_rt_1_1(all(~isnan(tr_rt_1_1),2),:); % remove NaNs form the vector
tr_ad_1_1 = tr_ad_1_1(all(~isnan(tr_ad_1_1),2),:);% remove NaNs form the vector
Vtr_rt_1_1 = Vtr_rt_1_1(all(~isnan(Vtr_rt_1_1),2),:);% remove NaNs form the vector
Vtr_ad_1_1 = Vtr_ad_1_1(all(~isnan(Vtr_ad_1_1),2),:);% remove NaNs form the vector
t_rt_1_1 = t_rt_1_1(all(~isnan(t_rt_1_1),2),:);% remove NaNs form the vector
t_ad_1_1 = t_ad_1_1(all(~isnan(t_ad_1_1),2),:);% remove NaNs form the vector

mean_tr_1_1 = 1/(t_rt_1_1(end)-t_rt_1_1(1)) .* trapz(t_rt_1_1, tr_rt_1_1);
if (numel(tr_ad_1_1>=0))
mean_ad_1_1 = 1/(t_ad_1_1(end)-t_ad_1_1(1)) .* trapz(t_ad_1_1, tr_ad_1_1);
end
