function calculate_abs_mag(cat,zl)
%%note: NO k correction applied this should be taken care of when doing real calculations!! 
file_data=load(cat);
[dlum,bla1,bla2]=lum_dist(zl);%in parsecs
abs_mag(:,1)=file_data(:,6)-5*(log10(dlum)-1);
abs_mag
counti=0;
noinbin(1:(25.5-14.5)/0.2+1)=0;
for i=-25.5:0.2:-14.5
    counti=counti+1;
       for t=1:length(file_data(:,1))
        if abs_mag(t)>i-0.5 && abs_mag(t)<i+0.5
        noinbin(counti)=noinbin(counti)+1;
        end
       end
    
end
hist(abs_mag)
median(abs_mag)
scatter(-25.5:0.2:-14.5,log(noinbin))
set(gca,'XDir','reverse')
file_data(:,12)=file_data(:,7);
file_data(:,7)=abs_mag(:,1);
save 'model_colors_allWAbsLum.txt' file_data -ascii
% %result is B band Vega
% file_data=load('model_colors_temp.txt');
% roundedup_z=round(zl*100);
% roundedup_z
% txitxo=load('SDSS_LRGv1.txt');
% for i=1:length(file_data(:,1))
% for t=1:length(txitxo(:,1))
%     if (txitxo(t,1)==roundedup_z/100)
%         L_0_solarBCG=txitxo(t,7)
%         m_0BCG=txitxo(t,4)
%     end
% end
% L_solarGal(i)=10^(L_0_solarBCG)*10^(-0.4*(file_data(i,6)-m_0BCG));
% end
% L_solarGal