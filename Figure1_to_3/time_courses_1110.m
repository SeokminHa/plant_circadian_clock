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

para=readmatrix('para_1110_full.csv');
binit=10^binding;

% toc1 time-series
ttoc1=[];stoc1=[];toc1ztld=[];toc1ztll=[];
toc1gi=[];toc1prr3=[];

% ztl time-series
tztl=[];sztld=[];sztll=[];ztldtoc1=[];ztlltoc1=[];
ztldgi=[];ztllgi=[];ztldprr3=[];ztllprr3=[];

%gi time-series
tgi=[];sgi=[];gitoc1=[];giztld=[];
giztll=[];giprr3=[];

% prr3 time-series
tprr3=[];sprr3=[];prr3toc1=[];prr3ztld=[];
prr3ztll=[];prr3gi=[];

% max values
max_toc1=[]; max_ztl=[]; max_gi=[]; max_prr3=[];
ztl_dark=[]; ztl_light=[];


size(para,1)
rel_amp=[];
toc1_stab=[];
toc1_stab_amp=[];

gp_p=[];
gp_g=[];
for jj=1:size(para,1)
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

    itdata=plevel(:,1)+plevel(:,6)+plevel(:,7)+plevel(:,8)+plevel(:,9); % total TOC1
    izdata=plevel(:,2)+plevel(:,3)+plevel(:,6)+plevel(:,7)+plevel(:,10)+plevel(:,11)+plevel(:,12)+plevel(:,13);  % total ZTL
    igdata=plevel(:,4)+plevel(:,8)+plevel(:,10)+plevel(:,11)+plevel(:,14); % total GI
    ip3data=plevel(:,5)+plevel(:,9)+plevel(:,12)+plevel(:,13)+plevel(:,14); % total PRR3

    tmax=max(itdata); zmax=max(izdata);
    gmax=max(igdata); pmax=max(ip3data);

    % normalize using maximum value
    itdata=itdata./tmax;  izdata=izdata./zmax;
    igdata=igdata./gmax;  ip3data=ip3data./pmax;

    % toc1
    ttoc1r=(plevel(:,1)+plevel(:,6)+plevel(:,7)+plevel(:,8)+plevel(:,9));
    stoc1r=plevel(:,1);
    toc1ztldr=plevel(:,6);
    toc1ztllr=plevel(:,7);
    toc1gir=plevel(:,8);
    toc1prr3r=plevel(:,9);

    % calculate the stability of TOC1
    toc1_rhythmic_degra = dt * (stoc1r+dtz1_z1 * toc1ztldr + dtz2_z2 * toc1ztllr + dtp_p * toc1prr3r)./ttoc1r;
    toc1_stability= 1 - toc1_rhythmic_degra/max(toc1_rhythmic_degra);
    toc1_stability_amp = max(toc1_stability)-min(toc1_stability);
    toc1_stab = [toc1_stab; toc1_stability'];
    toc1_stab_amp = [toc1_stab_amp;toc1_stability_amp];

    ttoc1r=ttoc1r/tmax';
    stoc1r=stoc1r./tmax';
    toc1ztldr=toc1ztldr./tmax';
    toc1ztllr=toc1ztllr./tmax';
    toc1gir=toc1gir./tmax';
    toc1prr3r=toc1prr3r./tmax';

    % ztl
    tztlr=plevel(:,2)+plevel(:,3)+plevel(:,6)+plevel(:,7)+plevel(:,10)+plevel(:,11)+plevel(:,12)+plevel(:,13);
    sztldr=plevel(:,2);
    sztllr=plevel(:,3);
    ztldtoc1r=plevel(:,6);
    ztlltoc1r=plevel(:,7);
    ztldgir=plevel(:,10);
    ztllgir=plevel(:,11);
    ztldprr3r=plevel(:,12);
    ztllprr3r=plevel(:,13);
    ztl_l = plevel(:,3)+plevel(:,7)+plevel(:,11)+plevel(:,13);
    ztl_d = plevel(:,2)+plevel(:,6)+plevel(:,10)+plevel(:,12);

    tztlr=tztlr./zmax;
    sztldr=sztldr./zmax;
    sztllr=sztllr./zmax;
    ztldtoc1r=ztldtoc1r./zmax;
    ztlltoc1r=ztlltoc1r./zmax;
    ztldgir=ztldgir./zmax;
    ztllgir=ztllgir./zmax;
    ztldprr3r=ztldprr3r./zmax;
    ztllprr3r=ztllprr3r./zmax;
    ztl_l = ztl_l./zmax;
    ztl_d = ztl_d./zmax;

    % gi
    tgir=plevel(:,4)+plevel(:,8)+plevel(:,10)+plevel(:,11)+plevel(:,14);
    sgir=plevel(:,4);
    gitoc1r=plevel(:,8);
    giztldr=plevel(:,10);
    giztllr=plevel(:,11);
    giprr3r=plevel(:,14);

    tgir=tgir./gmax;
    sgir=sgir./gmax;
    gitoc1r=gitoc1r./gmax;
    giztldr=giztldr./gmax;
    giztllr=giztllr./gmax;
    giprr3r=giprr3r./gmax;

    % prr3
    tprr3r=plevel(:,5)+plevel(:,9)+plevel(:,12)+plevel(:,13)+plevel(:,14);
    sprr3r=plevel(:,5);
    prr3toc1r=plevel(:,9);
    prr3ztldr=plevel(:,12);
    prr3ztllr=plevel(:,13);
    prr3gir=plevel(:,14);

    tprr3r=tprr3r./pmax;
    sprr3r=sprr3r./pmax;
    prr3toc1r=prr3toc1r./pmax;
    prr3ztldr=prr3ztldr./pmax;
    prr3ztllr=prr3ztllr./pmax;
    prr3gir=prr3gir./pmax;

    ttoc1=[ttoc1;ttoc1r'];stoc1=[stoc1;stoc1r'];
    toc1ztld=[toc1ztld;toc1ztldr'];toc1ztll=[toc1ztll;toc1ztllr'];
    toc1gi=[toc1gi;toc1gir'];toc1prr3=[toc1prr3;toc1prr3r'];
    max_toc1 = [max_toc1, tmax];

    tztl=[tztl;tztlr'];sztld=[sztld;sztldr'];sztll=[sztll;sztllr'];
    ztldtoc1=[ztldtoc1;ztldtoc1r'];ztlltoc1=[ztlltoc1;ztlltoc1r'];
    ztldgi=[ztldgi;ztldgir'];ztllgi=[ztllgi;ztllgir'];ztldprr3=[ztldprr3;ztldprr3r'];ztllprr3=[ztllprr3;ztllprr3r'];
    max_ztl = [max_ztl, zmax];
    ztl_light =[ztl_light; ztl_l'];  ztl_dark =[ztl_dark; ztl_d'];
    tgi=[tgi;tgir'];sgi=[sgi;sgir'];gitoc1=[gitoc1;gitoc1r'];giztld=[giztld;giztldr'];
    giztll=[giztll;giztllr'];giprr3=[giprr3;giprr3r']; gp_g=[gp_g,mean(giprr3r)];
    max_gi = [max_gi, gmax];

    tprr3=[tprr3; tprr3r'];sprr3=[sprr3; sprr3r'];prr3toc1=[prr3toc1;prr3toc1r'];prr3ztld=[prr3ztld;prr3ztldr'];
    prr3ztll=[prr3ztll;prr3ztllr'];prr3gi=[prr3gi;prr3gir']; gp_p=[gp_p,mean(prr3gir)];
    max_prr3 = [max_prr3, pmax];

    rel_amp= [rel_amp; [1-min(ttoc1r) 1-min(tztlr) 1-min(tgir) 1-min(tprr3r) ]];


end

csvwrite('TOC1 stability_amplitude.csv', toc1_stab_amp);


size(toc1_stability)

xlimit=[10^-3 10^3];
marsize=20;

tspan=0:0.1:24;


figure(11) % TOC1
subplot(2,3,1)
plot(tspan,mean(ttoc1),'r-','LineWidth', 2)
hold on

plot(tspan,mean(ttoc1)+std(ttoc1),'b-','LineWidth', 1)
plot(tspan,mean(ttoc1)-std(ttoc1),'r-','LineWidth', 1)
xlim([0 24])
ylim([0 1])
title('Total Graph for TOC1')
hold off

subplot(2,3,2)
plot(tspan,mean(toc1ztld),'m-','LineWidth', 2)
hold on
plot(tspan,mean(toc1ztld)+std(toc1ztld),'m-','LineWidth', 1)
plot(tspan,mean(toc1ztld)-std(toc1ztld),'m-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - ZTL_{dark} complex')

subplot(2,3,3)
plot(tspan,mean(toc1ztll),'m-','LineWidth', 2)
hold on
plot(tspan,mean(toc1ztll)+std(toc1ztll),'m-','LineWidth', 1)
plot(tspan,mean(toc1ztll)-std(toc1ztll),'m-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - ZTL_{light} complex')

subplot(2,3,4)
plot(tspan,mean(stoc1),'g-','LineWidth', 2)
hold on
plot(tspan,mean(stoc1)+std(stoc1),'g-','LineWidth', 1)
plot(tspan,mean(stoc1)-std(stoc1),'g-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('Single component TOC1')

subplot(2,3,5)
plot(tspan,mean(toc1gi),'y-','LineWidth', 2)
hold on
plot(tspan,mean(toc1gi)+std(toc1gi),'y-','LineWidth', 1)
plot(tspan,mean(toc1gi)-std(toc1gi),'y-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - GI complex')

subplot(2,3,6)
plot(tspan,mean(toc1prr3),'c-','LineWidth', 2)
hold on
plot(tspan,mean(toc1prr3)+std(toc1prr3),'c-','LineWidth', 1)
plot(tspan,mean(toc1prr3)-std(toc1prr3),'c-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - PRR3 complex')

figure(12) % ztl
subplot(3,4,1)
plot(tspan,mean(tztl),'r-','LineWidth', 1)
hold on
plot(tspan,mean(tztl)+std(tztl),'r-','LineWidth', 1)
plot(tspan,mean(tztl)-std(tztl),'r-','LineWidth', 1)
xlim([0 24])
ylim([0 1])
title('Total Graph for ZTL')
hold off

subplot(3,4,2)
plot(tspan,mean(ztl_light),'b-','LineWidth', 2,'DisplayName','ZTL_{light}')
hold on
plot(tspan,mean(ztl_dark),'r-','LineWidth', 2,'DisplayName','ZTL_{dark}')
hold on
plot(tspan,mean(ztl_dark)+std(ztl_dark),'r-','LineWidth', 1)
plot(tspan,mean(ztl_dark)-std(ztl_dark),'r-','LineWidth', 1)
plot(tspan,mean(ztl_light)+std(ztl_light),'k-','LineWidth', 1)
plot(tspan,mean(ztl_light)+std(ztl_light),'k-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('Total ZTL_{dark} and ZTL_{light}')

subplot(3,4,3)
plot(tspan,mean(sztld),'g-','LineWidth', 2)
hold on
plot(tspan,mean(sztld)+std(sztld),'g-','LineWidth', 1)
plot(tspan,mean(sztld)-std(sztld),'g-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('Single component ZTL_{dark}')

subplot(3,4,4)
plot(tspan,mean(sztll),'g-','LineWidth', 2)
hold on
plot(tspan,mean(sztll)+std(sztll),'g-','LineWidth', 1)
plot(tspan,mean(sztll)-std(sztll),'g-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('Single component ZTL_{light}')

subplot(3,4,5)
plot(tspan,mean(ztldtoc1),'m-','LineWidth', 2)
hold on
plot(tspan,mean(ztldtoc1)+std(ztldtoc1),'m-','LineWidth', 1)
plot(tspan,mean(ztldtoc1)-std(ztldtoc1),'m-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - ZTL_{dark} complex')

subplot(3,4,6)
plot(tspan,mean(ztlltoc1),'m-','LineWidth', 2)
hold on
plot(tspan,mean(ztlltoc1)+std(ztlltoc1),'m-','LineWidth', 1)
plot(tspan,mean(ztlltoc1)-std(ztlltoc1),'m-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('TOC1 - ZTL_{light} complex')

subplot(3,4,7)
plot(tspan,mean(ztldgi),'k-','LineWidth', 2)
hold on
plot(tspan,mean(ztldgi)+std(ztldgi),'k-','LineWidth', 1)
plot(tspan,mean(ztldgi)-std(ztldgi),'k-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title(' ZTL_{dark} - GI complex')

subplot(3,4,8)
plot(tspan,mean(ztllgi),'k-','LineWidth', 2)
hold on
plot(tspan,mean(ztllgi)+std(ztllgi),'k-','LineWidth', 1)
plot(tspan,mean(ztllgi)-std(ztllgi),'k-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title(' ZTL_{light} - GI complex')

subplot(3,4,9)
plot(tspan,mean(ztldprr3),'c-','LineWidth', 2)
hold on
plot(tspan,mean(ztldprr3)+std(ztldprr3),'c-','LineWidth', 1)
plot(tspan,mean(ztldprr3)-std(ztldprr3),'c-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('ZTL_{dark} - PRR3 complex')

subplot(3,4,10)
plot(tspan,mean(ztllprr3),'c-','LineWidth', 2)
hold on
plot(tspan,mean(ztllprr3)+std(ztllprr3),'c-','LineWidth', 1)
plot(tspan,mean(ztllprr3)-std(ztllprr3),'c-','LineWidth', 1)
hold off
xlim([0 24])
ylim([0 1])
title('ZTL_{light} - PRR3 complex')

figure(13) % GI
subplot(2,3,1)
plot(tspan,mean(tgi),'r-','LineWidth', 2)
hold on
plot(tspan,mean(tgi)+std(tgi),'r-','LineWidth', 1)
plot(tspan,mean(tgi)-std(tgi),'r-','LineWidth', 1)
xlim([0 24])
ylim([0 1])
title('Total Graph for GI')
hold off

subplot(2,3,2)
plot(tspan,mean(sgi),'g-','LineWidth', 2)
hold on;
plot(tspan,mean(sgi)+std(sgi),'g-','LineWidth', 1)
plot(tspan,mean(sgi)-std(sgi),'g-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('Single component GI')

subplot(2,3,3)
plot(tspan,mean(gitoc1),'y-','LineWidth', 2)
hold on;
plot(tspan,mean(gitoc1)+std(gitoc1),'y-','LineWidth', 1)
plot(tspan,mean(gitoc1)-std(gitoc1),'y-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('TOC1 - GI complex')

subplot(2,3,4)
plot(tspan,mean(giztld),'k-','LineWidth', 2)
hold on;
plot(tspan,mean(giztld)+std(giztld),'k-','LineWidth', 1)
plot(tspan,mean(giztld)-std(giztld),'k-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('ZTL_{dark} - GI complex')

subplot(2,3,5)
plot(tspan,mean(giztll),'k-','LineWidth', 2)
hold on;
plot(tspan,mean(giztll)+std(giztll),'k-','LineWidth', 1)
plot(tspan,mean(giztll)-std(giztll),'k-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('ZTL_{light} - GI complex')

subplot(2,3,6)
plot(tspan,mean(giprr3),'m-','LineWidth', 2)
hold on;
plot(tspan,mean(giprr3)+std(giprr3),'m-','LineWidth', 1)
plot(tspan,mean(giprr3)-std(giprr3),'m-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('GI - PRR3 complex')

figure(14) % PRR3
subplot(2,3,1)
plot(tspan,mean(tprr3),'r-','LineWidth', 2)
hold on
plot(tspan,mean(tprr3)+std(tprr3),'r-','LineWidth', 1)
plot(tspan,mean(tprr3)-std(tprr3),'r-','LineWidth', 1)
xlim([0 24])
ylim([0 1])
title('Total Graph for PRR3')
hold off

subplot(2,3,2)
plot(tspan,mean(sprr3),'g-','LineWidth', 2)
hold on;
plot(tspan,mean(sprr3)+std(sprr3),'g-','LineWidth', 1)
plot(tspan,mean(sprr3)-std(sprr3),'g-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('Single component PRR3')

subplot(2,3,3)
plot(tspan,mean(prr3toc1),'c-','LineWidth', 2)
hold on;
plot(tspan,mean(prr3toc1)+std(prr3toc1),'c-','LineWidth', 1)
plot(tspan,mean(prr3toc1)-std(prr3toc1),'c-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('TOC1 - PRR3 complex')

subplot(2,3,4)
plot(tspan,mean(prr3ztld),'k-','LineWidth', 2)
hold on;
plot(tspan,mean(prr3ztld)+std(prr3ztld),'k-','LineWidth', 1)
plot(tspan,mean(prr3ztld)-std(prr3ztld),'k-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('ZTL_{dark} - PRR3 complex')

subplot(2,3,5)
plot(tspan,mean(prr3ztll),'y-','LineWidth', 2)
hold on;
plot(tspan,mean(prr3ztll)+std(prr3ztll),'y-','LineWidth', 1)
plot(tspan,mean(prr3ztll)-std(prr3ztll),'y-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('ZTL_{light} - PRR3 complex')

subplot(2,3,6)
plot(tspan,mean(prr3gi),'m-','LineWidth', 2)
hold on;
plot(tspan,mean(prr3gi)+std(prr3gi),'m-','LineWidth', 1)
plot(tspan,mean(prr3gi)-std(prr3gi),'m-','LineWidth', 1)
hold off;
xlim([0 24])
ylim([0 1])
title('GI - PRR3 complex')

end


