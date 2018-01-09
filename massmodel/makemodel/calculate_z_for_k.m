function calculate_z_for_k(cluster_redshift, main_source_z, k_system);

% this function must have ad_dist.m and chiint.m in the same directory!!!

zc=cluster_redshift;



Dls1=ad_dist([zc main_source_z]);
Ds1=ad_dist(main_source_z);
% rat1=Dls(i)/Ds(i);
for i=1:500
second_source_z=0.4+(i-1)*0.01;   
Dls2=ad_dist([zc second_source_z]);
Ds2=ad_dist(second_source_z);
k_times(i)=(Dls2/Ds2)/(Dls1/Ds1);
end
k_dif=k_times-k_system;
k_dif=abs(k_dif);
min_dif=min(k_dif);
k_ind=find(k_dif==min_dif);
z_right=0.4+(k_ind-1)*0.01