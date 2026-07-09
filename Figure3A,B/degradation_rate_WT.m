function time_courses_1110


global tt tz tg tp ...
    dt dz dg dp ...
    dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p ...
    kc1 kc2 ... %AB -> A : dAB_A*dB
    bb ...
    ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp ...
    light  days toc1mrna gimrna prr3mrna

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

binding=1;

bd_type=[1,1,1,0];
itdata_tot=[];
toc1_rhythmic_degra_tot=[];
para=readmatrix('para_1110.csv');
binit=10^binding;


size(para,1)
toc1_stab=[];
toc1_stab_amp=[];

gp_p=[];
gp_g=[];
for jj=1:10%size(para,1)
    if rem(jj,20)==0
        jj
    end
    paraset=para(jj,:);
    tranp=paraset(1:4);
    degp=paraset(5:8);
    conp=paraset(9:10);
    unbindp=paraset(11:19);
    ratep=paraset(20:37);

    tmp=num2cell(tranp); dmp=num2cell(degp); kmp=num2cell(conp);
    ubmp=num2cell(unbindp); cdmp=num2cell(ratep);

    [tt tz tg tp]=deal(tmp{:});
    [dt dz dg dp]=deal(dmp{:});
    [kc1 kc2]=deal(kmp{:});
    [bb]=binit;
    [ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp]=deal(ubmp{:});
    [dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p]=deal(cdmp{:});
    days=7;
    plevel = [];
    C1=0*ones(1,14);
    for j=1:days % Simulate using ODE
        light = 1; %ZT0~ZT12, light
        tspan = 24*(j-1):0.1:24*(j-1)+12;
        options = odeset('NonNegative',1:14);
        [T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
        if j==days
            plevel = [plevel; C1];
        end

        light = 0; %ZT12~ZT24, no light
        tspan = 24*(j-1)+12:0.1:24*j;
        options = odeset('NonNegative',1:14);
        [T1,C1] = ode15s(@(t,C) multi_degradation_ODE_v6_5(t,C,bd_type),tspan,C1(end,:),options);
        if j==days
            plevel=[plevel; C1(2:end,:)];
        end
    end
        
    % toc1
    ttoc1r=(plevel(:,1)+plevel(:,6)+plevel(:,7)+plevel(:,8)+plevel(:,9));
    stoc1r=plevel(:,1);
    toc1ztldr=plevel(:,6);
    toc1ztllr=plevel(:,7);
    toc1gir=plevel(:,8);
    toc1prr3r=plevel(:,9);

    % calculate the degradation rate of TOC1
    toc1_rhythmic_degra = dt * (stoc1r+dtz1_z1 * toc1ztldr + dtz2_z2 * toc1ztllr + dtp_p * toc1prr3r)./ttoc1r;
    toc1_rhythmic_degra_tot = [toc1_rhythmic_degra_tot; toc1_rhythmic_degra'];
end

% csvwrite('TOC1 degra_wt.csv', toc1_rhythmic_degra_tot);
end


