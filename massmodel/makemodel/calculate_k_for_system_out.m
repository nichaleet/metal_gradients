function [dlsds_ratio]=calculate_k_for_system_out(cluster_redshift, main_source_z, second_source_z);

% this function must have ad_dist.m and chiint.m in the same directory!!!

zc=cluster_redshift;



Dls1=ad_dist([zc main_source_z]);
Ds1=ad_dist(main_source_z);
% rat1=Dls(i)/Ds(i);
Dls2=ad_dist([zc second_source_z]);
Ds2=ad_dist(second_source_z);
k_times=(Dls2/Ds2)/(Dls1/Ds1)
dlsds_ratio=k_times;