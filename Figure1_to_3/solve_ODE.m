function [imax, zd_in_l, zl_in_d, itd, izd, igd, ip3d, stab_less, stab_more, zt1, zt13] = solve_ODE(parav1, parav2, parav3, parav4, parav5, bd_type, binit)

global tt tz tg tp ...
    dt dz dg dp ...
    dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p ... %AB -> A : dAB_A*dB
    kc1 kc2 ...
    bb ...
    ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp ...
    light  days toc1mrna gimrna prr3mrna...
    tbestfit zbestfit gbestfit pbestfit maxfit...
    value temp iteration
toc1mrna=[0 1 5 9 13 17 21 24; ...
    0.401508 0.376 0.376 0.69 1 0.52 0.489 0.401508];
gimrna=[0 3 6 9 12 15 18 21 24; ...
    0.0535789, 0.277942, 0.813305, 1., 0.373043, 0.00648925, 0.00439222, 0.0122333, 0.0535789];
prr3mrna=[0 3 6 9 12 15 18 21 24; ...
    0.010205, 0.00916596, 0.126271, 0.801952, 1., 0.091304, 0.0357569, 0.022007, 0.010205];

toc1p=[1 5 9 13 17 21; ...
    0.0649 0.0346 0.29 0.987 1 0.645];
ztlp=[1, 5, 9, 13, 17, 21; ...
    0.115, 0.187, 0.445, 1., 0.718, 0.56];
gip=[0 3 6 9 12 15 18 21 24; ...
    0.237939, 0.0842713, 0.365812, 0.913379, 1., 0.425148, 0.208709, 0.0937085, 0.096325];
prr3p=[0 3 6 9 12 15 18 21 24; ...
    0.021049, 0.0711328, 0.128753, 0.574524, 1., 0.587505, 0.371859, 0.355726, 0.104436];

tmp=num2cell(parav1); dmp=num2cell(parav2); kmp=num2cell(parav3);
ubmp=num2cell(parav4); cdmp=num2cell(parav5);

[tt tz tg tp]=deal(tmp{:});
[dt dz dg dp]=deal(dmp{:});
[kc1 kc2]=deal(kmp{:});
[bb]=binit;
[ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp]=deal(ubmp{:});
[dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p]=deal(cdmp{:});

days=7;
plevel = [];
C1=0*ones(1,14);

for j=1:days
    light = 1;
    tspan = 24*(j-1):1:24*(j-1)+12;
    options = odeset('NonNegative',1:14);
    [T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
    if j==days
        plevel = [plevel; C1];
    end
    
    light = 0;
    tspan = 24*(j-1)+12:1:24*j;
    options = odeset('NonNegative',1:14);
    [T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
    if j==days
        plevel=[plevel; C1(2:end,:)];
    end
end

fin=C1(end,:);
tt1=tt;
tz1=tz;
tg1=tg;
tp1=tp;
light = 1;
tspan = 24*days:0.5:24*days+1;
options = odeset('NonNegative',1:14);
[T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,fin,options);
z0=C1(end,2)+C1(end,6)+C1(end,10)+C1(end,12)+C1(end,3)+C1(end,7)+C1(end,11)+C1(end,13);
tt=0;
tz=0;
tg=0;
tp=0;
tspan = 24*days+1:0.5:24*days+2;
options = odeset('NonNegative',1:14);
[T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
z1tot=C1(end-1,:);
z2tot=C1(end,:);

z1=z1tot(2)+z1tot(6)+z1tot(10)+z1tot(12)+z1tot(3)+z1tot(7)+z1tot(11)+z1tot(13);
z2=z2tot(2)+z2tot(6)+z2tot(10)+z2tot(12)+z2tot(3)+z2tot(7)+z2tot(11)+z2tot(13);

zt1 = [z1/z0,z2/z0];

tt=tt1;
tz=tz1;
tg=tg1;
tp=tp1;

light = 1;
tspan = 24*days:0.5:24*days+12;
options = odeset('NonNegative',1:14);
[T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,fin,options);
light=0;
tspan = 24*days+12:0.5:24*days+13;
options = odeset('NonNegative',1:14);
[T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
zz0=C1(end,2)+C1(end,6)+C1(end,10)+C1(end,12)+C1(end,3)+C1(end,7)+C1(end,11)+C1(end,13);
tt=0;
tz=0;
tg=0;
tp=0;
tspan = 24*days+13:0.5:24*days+14;
options = odeset('NonNegative',1:14);
[T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
z1tot=C1(end-1,:);
z2tot=C1(end,:);
zz1=z1tot(2)+z1tot(6)+z1tot(10)+z1tot(12)+z1tot(3)+z1tot(7)+z1tot(11)+z1tot(13);
zz2=z2tot(2)+z2tot(6)+z2tot(10)+z2tot(12)+z2tot(3)+z2tot(7)+z2tot(11)+z2tot(13);

zt13 = [zz1/zz0,zz2/zz0];
tt=tt1;
tz=tz1;
tg=tg1;
tp=tp1;

itdata=plevel(:,1)+plevel(:,6)+plevel(:,7)+plevel(:,8)+plevel(:,9);
izdata=plevel(:,2)+plevel(:,3)+plevel(:,6)+plevel(:,7)+plevel(:,10)+plevel(:,11)+plevel(:,12)+plevel(:,13);
igdata=plevel(:,4)+plevel(:,8)+plevel(:,10)+plevel(:,11)+plevel(:,14);
ip3data=plevel(:,5)+plevel(:,9)+plevel(:,12)+plevel(:,13)+plevel(:,14);

iz1data=plevel(:,2)+plevel(:,6)+plevel(:,10)+plevel(:,12);
iz2data=plevel(:,3)+plevel(:,7)+plevel(:,11)+plevel(:,13);
zd_in_l =  mean(iz1data(1:13)) / mean(izdata(1:13));
zl_in_d = mean(iz2data(13:25)) / mean(izdata(13:25));

itdata=itdata(toc1p(1,:)+1,1);  izdata=izdata(ztlp(1,:)+1,1);
igdata=igdata(gip(1,:)+1,1);  ip3data=ip3data(prr3p(1,:)+1,1);

tmax = max(itdata); zmax = max(izdata);
gmax = max(igdata); pmax = max(ip3data);
imax = [tmax, zmax, gmax, pmax];

itd=itdata./tmax;  izd=izdata./zmax;
igd=igdata./gmax;  ip3d=ip3data./pmax;
% 
tzd_amount1 = plevel(:,6)/tmax;
tzd_amount1 = tzd_amount1(13:25);
tzd_night1 = mean( tzd_amount1 );
% 
% tzd_amount2 = plevel(:,6)/zmax;
% tzd_amount2 = tzd_amount2(13:25);
% tzd_night2 = mean( tzd_amount2 );

tzl_amount1 = plevel(:,7)/tmax;
tzl_amount1 = tzl_amount1(1:13);
tzl_day1 = mean( tzl_amount1 );

% tzl_amount2 = plevel(:,7)/zmax;
% tzl_amount2 = tzl_amount2(1:13);
% tzl_day2 = mean( tzl_amount2 );

tp_amount1 = plevel(:,9)/tmax;
tp1 = mean( tp_amount1 );

% tp_amount2 = plevel(:,9)/pmax;
% tp2 = mean( tp_amount2 );

zlg_amount1 = plevel(:,11)/zmax;
zlg_amount1 = zlg_amount1(1:13);
zlg_day1 = mean( zlg_amount1 );

zlg_amount2 = plevel(:,11)/gmax;
zlg_amount2 = zlg_amount2(1:13);
zlg_day2 = mean( zlg_amount2 );

stab_less = [tzd_night1,  tzl_day1];
% stab_less = [tzd_night1, tzd_night2];
stab_more = [tp1, zlg_day1, zlg_day2 ];




end