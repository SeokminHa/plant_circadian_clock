tt_degra_wt1 = readmatrix('TOC1 degra_wt.csv');
tt_degra_z1 = readmatrix('TOC1 degra_z.csv');
tt_degra_p1 = readmatrix('TOC1 degra_p.csv');
tt_degra_zp1 = readmatrix('TOC1 degra_zp.csv');

time_wt =tt_degra_wt1(1,:);
tt_degra_wt = tt_degra_wt1(2:end,:);

time_z =tt_degra_z1(1,:);
tt_degra_z = tt_degra_z1(2:end,:);

time_p =tt_degra_p1(1,:);
tt_degra_p = tt_degra_p1(2:end,:);

time_zp =tt_degra_zp1(1,:);
tt_degra_zp = tt_degra_zp1(2:end,:);


% Figure 3A, stability time course
csvwrite('TOC1_raw_stab_wt_summary.csv', [time_wt;3-mean(tt_degra_wt);std(tt_degra_wt)]);
csvwrite('TOC1_raw_stab_z_summary.csv', [time_z;3-mean(tt_degra_z);std(tt_degra_z)]);
csvwrite('TOC1_raw_stab_p_summary.csv', [time_p;3-mean(tt_degra_p);std(tt_degra_p)]);
csvwrite('TOC1_raw_stab_zp_summary.csv', [time_zp;3-mean(tt_degra_zp);std(tt_degra_zp)]);


% Figure 3B, amplitude   
stab_wt_mean0 = max(tt_degra_wt, [], 2) - min(tt_degra_wt, [], 2);
stab_z_mean0  = max(tt_degra_z, [], 2)  - min(tt_degra_z, [], 2);
stab_p_mean0  = max(tt_degra_p, [], 2)  - min(tt_degra_p, [], 2);
stab_zp_mean0 = max(tt_degra_zp, [], 2) - min(tt_degra_zp, [], 2);

csvwrite('TOC1_stab_amp_raw_all.csv', [ stab_wt_mean0, stab_p_mean0, stab_z_mean0 ,stab_zp_mean0 ]);
