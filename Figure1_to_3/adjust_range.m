function [parav1, parav2, parav3, parav4, parav5] = adjust_range(bd_type, transmax, transmin, degra, parav1, parav2, parav3, parav4, parav5, th_par, th1, th2, m)

global tt tz tg tp ...
    dt dz dg dp ...
    dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p ... %AB -> A : dAB_A*dB
    kc1 kc2 ...
    bb ...
    ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp ...
    light  days toc1mrna gimrna prr3mrna...
    tbestfit zbestfit gbestfit pbestfit maxfit...
    value temp iteration

for i = [1 2 4 6] % bind well
    parav4(i) = max(10^(m-3), parav4(i));
    parav4(i) = min(parav4(i), 10^(m-3+th1));
end
for i = [5] % not bind well
    parav4(i) = max(10^(m+3-th1), parav4(i));
    parav4(i) = min(parav4(i), 10^(m+3));
end

parav5(parav5>10^1.31)=10^1.31;
parav5(parav5 > 0 & parav5<10^-1.31)=10^-1.31;
for j= [2 4] % de-stabilize
    parav5(j) = max(10^th2, parav5(j));
end
for i= [1 3 7] % no effect
    parav5(i) = 1;% max(10^(th_par/log10(parav5(i+1))), parav5(i) );
end
for j= [8 13 14] % stabilize
    parav5(j) = max(10^-1.31, parav5(j));
    parav5(j) = min(parav5(j), 10^-th2);
end
for j= [11 15 17] % stabilize / new finding force
    parav5(j) = max(10^-1.31, parav5(j));
    parav5(j) = min(parav5(j), 10^-0.1);
end
for j= [12 16 18] % destabilize / new finding force
    parav5(j) = max(10^0.1, parav5(j));
    parav5(j) = min(parav5(j), 10^1.31);
end

% GI and PRR3 bind
    parav4(9) = max(10^(m-3), parav4(9));
    parav4(9) = min(parav4(9), 10^(m));

% ZTL_dark and PRR bind
    parav4(7) = max(10^(m-3), parav4(7));
    parav4(7) = min(parav4(7), 10^(m-3));

% ZTL_light and PRR3 bind
    parav4(8) = max(10^(m-3), parav4(8));
    parav4(8) = min(parav4(8), 10^(m-3));

% TOC1 and GI not bind
    parav4(3) = max(10^(m), parav4(3));
    parav4(3) = min(parav4(3), 10^(m+3));

% for i=[17 11 15 5]
%     if parav5(i+1) >=1
%         parav5(i) = max(10^(max(-1.31, th_par/log10(parav5(i+1)))), parav5(i) );
%     else
%         parav5(i) = min(10^( min(1.31, th_par/log10(parav5(i+1)))), parav5(i) );
%     end
% end
if bd_type(1) ==0
    parav4(9) = 0;
    parav5(17) = 0;
    parav5(18) = 0;
end
if bd_type(2) ==0
    parav4(7) = 0;
    parav5(11) = 0;
    parav5(12) = 0;
end
if bd_type(3) ==0
    parav4(8) = 0;
    parav5(15) = 0;
    parav5(16) = 0;
end
if bd_type(4) ==0
    parav4(3) = 0;
    parav5(5) = 0;
    parav5(6) = 0;
end

parav1(parav1>transmax)=transmax;
parav1(parav1<transmin)=transmin;
parav2(parav2>degra)=degra;
parav2(parav2<0)=0;
parav3(parav3>100)=100;
parav3(parav3<0)=0;
parav4(parav4>10^(m+3))=10^(m+3);
parav4(parav4 > 0 & parav4<10^(m-3))=10^(m-3);
parav5(parav5>10^1.31)=10^1.31;
parav5(parav5 > 0 & parav5<10^-1.31)=10^-1.31;

end