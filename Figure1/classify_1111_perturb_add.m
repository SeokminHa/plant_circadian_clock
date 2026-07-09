function classify_1111_perturb_add


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


binding=1;

upper = 10^0.1;
lower = 10^-0.1;

error_thr = 0.5013;

cat_tot = [];
para=readmatrix('Supplementary Data 2_2000.csv');

binit=10^binding;


test_bd_types = [0, 1, 1, 1; ...
                 1, 0, 1, 1; ...
                 1, 1, 0, 1; ...
                 1, 1, 1, 0];

perturb_data = [];

for jj = 1:size(para, 1)
    fprintf('Iteration: %d / %d\n', jj, size(para, 1));

    paraset = para(jj, :);
    tranp = paraset(1:4);
    degp = paraset(5:8);
    conp = paraset(9:10);
    unbindp = paraset(11:19);
    ratep = paraset(20:37);

    row_errors = [];

    for k = 1:4
        current_bd_type = test_bd_types(k, :);

        [imax, zd_in_l, zl_in_d, itd, izd, igd, ip3d, stab_less, stab_more, zt1, zt13] = ...
            solve_ODE(tranp, degp, conp, unbindp, ratep, current_bd_type, binit);

        err = max([norm(toc1p(2,:)-itd'), norm(ztlp(2,:)-izd'), norm(gip(2,:)-igd'), norm(prr3p(2,:)-ip3d')]) ...
              + weight_chx1 * norm([m11 m12]-zt1) + weight_chx2 * norm([m21 m22]-zt13);
        err_t = norm(toc1p(2,:)-itd');
        err_z = norm(ztlp(2,:)-izd');
        err_g = norm(gip(2,:)-igd');
        err_p = norm(prr3p(2,:)-ip3d');
        err_chx1 = norm([m11 m12]-zt1);
        err_chx2 = norm([m21 m22]-zt13);

        row_errors = [row_errors, err, err_t, err_z, err_g, err_p, err_chx1, err_chx2];
    end

    perturb_data = [perturb_data; row_errors];

end


para_1110 = [];

for jj=1:size(para,1)

    cat = ones(1,8);

    paraset=para(jj,:);
    unbindp=paraset(11:19);
    ratep=paraset(20:37);

    % TOC1 - GI
    if unbindp(3) < 10 && perturb_data(jj, 4) > error_thr
        cat(1) = 4;
        cat(2) = 4;
        if ratep(5) > upper
            cat(1) = 3;
        elseif ratep(5) < lower
            cat(1) = 2;
        end
        if ratep(6) > upper
            cat(2) = 3;
        elseif ratep(6) < lower
            cat(2) = 2;
        end
    end
    % ZTL - PRR3
    if unbindp(7) < 10 && perturb_data(jj, 2) > error_thr
        cat(3) = 4;
        cat(4) = 4;
        if ratep(12) > upper
            cat(3) = 3;
        elseif ratep(12) < lower
            cat(3) = 2;
        end
        if ratep(11) > upper
            cat(4) = 3;
        elseif ratep(11) < lower
            cat(4) = 2;
        end
    end
    % ZTL - PRR3
    if unbindp(8) < 10 && perturb_data(jj, 3) > error_thr
        cat(5) = 4;
        cat(6) = 4;
        if ratep(15) > upper
            cat(5) = 3;
        elseif ratep(15) < lower
            cat(5) = 2;
        end
        if ratep(16) > upper
            cat(6) = 3;
        elseif ratep(16) < lower
            cat(6) = 2;
        end
    end

    % GI - PRR3
    if unbindp(9) < 10 && perturb_data(jj, 1) > error_thr
        cat(7) = 4;
        cat(8) = 4;
        if ratep(17) > upper
            cat(7) = 3;
        elseif ratep(17) < lower
            cat(7) = 2;
        end
        if ratep(18) > upper
            cat(8) = 3;
        elseif ratep(18) < lower
            cat(8) = 2;
        end
    end

    cat_tot =[cat_tot;cat];
    if sum(cat == [1,1,2,2,2,2,2,2]) == 8
        para_1110 = [para_1110; paraset];

    end
end
size(para_1110,1)
csvwrite('para_1110_full.csv', para_1110);

csvwrite('classify_result_1111_essential.csv', cat_tot);

[unique_rows, ~, idx] = unique(cat_tot, 'rows');
counts = accumarray(idx, 1);
[counts_sorted, order] = sort(counts, 'descend');
rows_sorted = unique_rows(order, :);
result = [rows_sorted, counts_sorted];
csvwrite('classify_result_sorted_1111_essential.csv', result);


end


