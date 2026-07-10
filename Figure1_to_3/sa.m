% Full model SA for bb = 10^1
% assume degradation rate of ZTL dark and light is same
% impossible for A -> B stabilize & B -> A destabilize
% ZTL dark in light, light in dark < 20%
% TOC1-ZTL complex less than <15% of total TOC1 in dark phase
% TOC1-ZTL complex less than <15% of total ZTL in dark phase

function sa(bd_type, ran_num)

rng(ran_num);

global tt tz tg tp ...
    dt dz dg dp ...
    dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p ... %AB -> A : dAB_A*dB
    kc1 kc2 ...
    bb ...
    ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp ...
    light  days toc1mrna gimrna prr3mrna...
    value temp iteration error error_t error_z error_p error_g  error_chx1 error_chx2

th_stab=10;
transmin=0.25;
transmax=4;

result1=[]; result2=[];

z1min1=0.3483;
z1max1=0.5824;
m11=(z1min1+z1max1)/2;

z1min2=0.2022;
z1max2=0.290;
m12=(z1min2+z1max2)/2;


z13min1=0.6737;
z13max1=0.8416;
m21=(z13min1+z13max1)/2;

z13min2=0.5634;
z13max2=0.7524;
m22=(z13min2+z13max2)/2;


binit= 10; % we fix binding rate to be 10

% mRNA data
toc1mrna=[0 1 5 9 13 17 21 24; ...
    0.401508 0.376 0.376 0.69 1 0.52 0.489 0.401508];
gimrna=[0 3 6 9 12 15 18 21 24; ...
    0.0535789, 0.277942, 0.813305, 1., 0.373043, 0.00648925, 0.00439222, 0.0122333, 0.0535789];
prr3mrna=[0 3 6 9 12 15 18 21 24; ...
    0.010205, 0.00916596, 0.126271, 0.801952, 1., 0.091304, 0.0357569, 0.022007, 0.010205];

% protein data
toc1p=[1 5 9 13 17 21; ...
    0.0649 0.0346 0.29 0.987 1 0.645];
ztlp=[1, 5, 9, 13, 17, 21; ...
    0.115, 0.187, 0.445, 1., 0.718, 0.56];
gip=[0 3 6 9 12 15 18 21 24; ...
    0.237939, 0.0842713, 0.365812, 0.913379, 1., 0.425148, 0.208709, 0.0937085, 0.096325];
prr3p=[0 3 6 9 12 15 18 21 24; ...
    0.021049, 0.0711328, 0.128753, 0.574524, 1., 0.587505, 0.371859, 0.355726, 0.104436];

degra = 2; % degradation rate: [0, degra]
val = 10;
th_par = -0.1;

weight_chx1 = 0.5; 
weight_chx2 = 1; 

for wb=1:200
    
    value=100000;
    testmind=8;
    testmind_un=18;
    iteration=0;
    temp=1;
    
    t5=0; m=1;
    
    itran=transmin+(transmax-transmin)*rand(1,4);
    idran=0+(degra-0)*rand(1,4);
    ikran=0+(100-0)*rand(1,2);
    iuran=10.^(m-3+ 6*rand(1,9));
    icran=10.^(-1.31+(1.31-(-1.31))*rand(1,18));
    for zz = [1 3 7]
        icran(zz)=1;
    end
    th1=1;
    th2=0.5;
    [itran, idran, ikran, iuran, icran] = adjust_range(bd_type, transmax, transmin, degra, itran, idran, ikran, iuran, icran, th_par, th1, th2, m);
    tir=itran; dir=idran; kir=ikran; ubir=iuran; cdir = icran;
    rr=[tir,dir,kir,binit,ubir,cdir];
    
    maxp=0.4;
    hold=0.35;
    
    while iteration < 10001 && value > hold
                if value < 0.4 & t5 == 0
                    [result1_5,result2_5] = save_to_csv(bd_type, degra, ran_num, 4, result1_5,result2_5);
                    t5=1;
                end
        
        thres1=max([-maxp -value/testmind]);
        thres2=min([maxp value/testmind]);
        
        uthres1=max([-maxp -value/testmind_un]);
        uthres2=min([maxp value/testmind_un]);
        
        perturb1=exp(thres1+(thres2-thres1)*rand(1,4));
        perturb2=exp(thres1+(thres2-thres1)*rand(1,4));
        perturb3=exp(thres1+(thres2-thres1)*rand(1,2));
        perturb4=exp(uthres1+(uthres2-uthres1)*rand(1,9));
        perturb5=exp(thres1+(thres2-thres1)*rand(1,18));
        
        parav1=perturb1.*tir; parav2=perturb2.*dir;
        parav3=perturb3.*kir; parav5=perturb5.*cdir;
        
        ublogscale=log10(ubir)+(3-m);
        parav4=10.^(perturb4.*ublogscale - (3-m));
        
        [parav1, parav2, parav3, parav4, parav5] = adjust_range(bd_type, transmax, transmin, degra, parav1, parav2, parav3, parav4, parav5, th_par, th1, th2, m);
        
        [imax, zd_in_l, zl_in_d, itd, izd, igd, ip3d, stab_less, stab_more, zt1, zt13] = solve_ODE(parav1, parav2, parav3, parav4, parav5, bd_type, binit);
        while ~(max(imax)/min(imax) < val &   zd_in_l < 0.15 & zl_in_d < 0.15  & ...
                size(stab_more,2) == size(stab_more(stab_more>th_stab/100),2) & size(stab_less,2) == size(stab_less(stab_less<th_stab/100),2) )
           
            iteration
            
            [max(imax)/min(imax) size(stab_more(stab_more>th_stab/100),2)+size(stab_less(stab_less<th_stab/100),2) max(zd_in_l,zl_in_d)]
            
            thres1=max([-maxp -value/testmind]);
            thres2=min([maxp value/testmind]);
            
            uthres1=max([-maxp -value/testmind_un]);
            uthres2=min([maxp value/testmind_un]);
            
            perturb1=exp(thres1+(thres2-thres1)*rand(1,4));
            perturb2=exp(thres1+(thres2-thres1)*rand(1,4));
            perturb3=exp(thres1+(thres2-thres1)*rand(1,2));
            perturb4=exp(uthres1+(uthres2-uthres1)*rand(1,9));
            perturb5=exp(thres1+(thres2-thres1)*rand(1,18));
            
            if iteration ==0
                if ~( max(imax)/min(imax) < val)
                    tir = transmin+(transmax-transmin)*rand(1,4);
                    dir = 0+(degra-0)*rand(1,4);
                end
                if ~(zd_in_l < 0.15 & zl_in_d < 0.15)
                    perturb3(1) = max([ 1/perturb3(1) perturb3(1)]);
                    perturb3(2) = max([ 1/perturb3(2) perturb3(2)]);
                end
                kir = kir.*perturb3;
                ubir=10.^(m-3+ 6*rand(1,9));
                cdir=10.^(-1.31+(1.31-(-1.31))*rand(1,18));
                
                [tir, dir, kir, ubir, cdir] = adjust_range(bd_type, transmax, transmin, degra,tir, dir, kir, ubir, cdir, th_par, th1, th2, m);
                kir
                parav1=tir; parav2=dir;
                parav3=kir; parav5=cdir;
                
                ublogscale=log10(ubir)+(3-m);
                parav4=10.^(perturb4.*ublogscale - (3-m));
                
                [parav1, parav2, parav3, parav4, parav5] = adjust_range(bd_type, transmax, transmin, degra, parav1, parav2, parav3, parav4, parav5, th_par, th1, th2, m);
                [imax, zd_in_l, zl_in_d, itd, izd, igd, ip3d, stab_less, stab_more, zt1, zt13] = solve_ODE(parav1, parav2, parav3, parav4, parav5, bd_type, binit);
            else % iteration > 0
                if ~( max(imax)/min(imax) < val)
                    max_index = find( imax == max(imax) );
                    min_index = find( imax == min(imax) );
                    perturb1(max_index) = min([ 1/perturb1(max_index) perturb1(max_index)]);
                    perturb1(min_index) = max([1/perturb1(min_index) perturb1(min_index)]);
                end
                parav1=perturb1.*tir;
                parav2=perturb2.*dir;
                parav3=perturb3.*kir;
                ublogscale=log10(ubir)+(3-m);
                parav4=10.^(perturb4.*ublogscale - (3-m));
                parav5=perturb5.*cdir;
                
                [parav1, parav2, parav3, parav4, parav5] = adjust_range(bd_type, transmax, transmin, degra, parav1, parav2, parav3, parav4, parav5, th_par, th1, th2, m);
                [imax, zd_in_l, zl_in_d, itd, izd, igd, ip3d,  stab_less, stab_more, zt1, zt13] = solve_ODE(parav1, parav2, parav3, parav4, parav5, bd_type, binit);
            end
            
        end
        error_pre = max( [norm(toc1p(2,:)-itd'), norm(ztlp(2,:)-izd'), norm(gip(2,:)-igd'),  norm(prr3p(2,:)-ip3d')] ) + weight_chx1 * norm([m11 m12]-zt1) + weight_chx2 * norm([m21 m22]-zt13);
        
        temp=0.999*temp;
        iteration=iteration+1;
        
        if error_pre < value
            error = max( [norm(toc1p(2,:)-itd'), norm(ztlp(2,:)-izd'), norm(gip(2,:)-igd'),  norm(prr3p(2,:)-ip3d')] ) + weight_chx1 * norm([m11 m12]-zt1) + weight_chx2 * norm([m21 m22]-zt13);
            error_t =  norm(toc1p(2,:)-itd');
            error_z = norm(ztlp(2,:)-izd');
            error_g = norm(gip(2,:)-igd');
            error_p = norm(prr3p(2,:)-ip3d');
            error_chx1 = norm([m11 m12]-zt1);
            error_chx2 = norm([m21 m22]-zt13);            
            value=error;
            paraset1=[tt tz tg tp];
            paraset2=[dt dz dg dp];
            paraset3=[kc1 kc2];
            paraset4=[ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp];
            paraset5=[dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p];
            
            tir=paraset1; dir=paraset2; kir=paraset3; ubir=paraset4; cdir=paraset5;
            
            disp([value, sqrt( error_t^2+ error_z^2+ error_g^2+ error_p^2 )  ,error_chx2, log10(dgp_g), iteration/10000])
            
        elseif exp(-abs(value-error_pre)/temp) > rand
            error = max( [norm(toc1p(2,:)-itd'), norm(ztlp(2,:)-izd'), norm(gip(2,:)-igd'),  norm(prr3p(2,:)-ip3d')] ) + weight_chx1 * norm([m11 m12]-zt1) + weight_chx2 * norm([m21 m22]-zt13);
            error_t =  norm(toc1p(2,:)-itd');
            error_z = norm(ztlp(2,:)-izd');
            error_g = norm(gip(2,:)-igd');
            error_p = norm(prr3p(2,:)-ip3d');
            error_chx1 = norm([m11 m12]-zt1);
            error_chx2 = norm([m21 m22]-zt13);            
            value=error;
            paraset1=[tt tz tg tp];
            paraset2=[dt dz dg dp];
            paraset3=[kc1 kc2];
            paraset4=[ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp];
            paraset5=[dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p];
            
            tir=paraset1; dir=paraset2; kir=paraset3; ubir=paraset4; cdir=paraset5;
            
            disp([1,value, max( [error_t, error_z, error_g, error_p]), error_chx2, iteration/10000])
        else
            disp([2, value, max( [error_t, error_z, error_g, error_p]) ,error_chx2])
        end
    end
    [result1,result2] = save_to_csv(bd_type, tir, dir, kir, ubir, cdir, degra, ran_num, 0, weight_chx1, weight_chx2, result1, result2);
    
end

end