function [result1,result2]=save_to_csv(bd_type, parav1, parav2, parav3, parav4, parav5, degra, ran_num, th, weight_chx1,  weight_chx2, result1,result2)
global tt tz tg tp ...
    dt dz dg dp ...
    dtz1_t dtz1_z1 dtz2_t dtz2_z2 dtg_t dtg_g dtp_t dtp_p dz1g_z1 dz1g_g dz1p_z1 dz1p_p dz2g_z2 dz2g_g dz2p_z2 dz2p_p dgp_g dgp_p ... %AB -> A : dAB_A*dB
    kc1 kc2 ...
    ubtz1 ubtz2 ubtg ubtp ubzg1 ubzg2 ubzp1 ubzp2 ubgp ...
    value temp iteration error_t error_z error_p error_g error_chx1 error_chx2

            result1=[result1; [parav1, parav2, parav3, parav4, parav5]];
            result2=[result2; [value, temp, iteration error_t error_z error_g error_p error_chx1 error_chx2]];
            csvwrite(['para_',num2str(bd_type(1)),num2str(bd_type(2)),num2str(bd_type(3)),num2str(bd_type(4)),'_',num2str(degra),'_',num2str(weight_chx1),'_',num2str(weight_chx2),'_',num2str(ran_num),'_0.',num2str(th),'.csv'],result1)
            csvwrite(['error_',num2str(bd_type(1)),num2str(bd_type(2)),num2str(bd_type(3)),num2str(bd_type(4)),'_',num2str(degra),'_',num2str(weight_chx1),'_',num2str(weight_chx2),'_',num2str(ran_num),'_0.',num2str(th),'.csv'],result2)  
            
end