function [magnification_ALL]=calculate_magnification_ALL_highz;
%This function calculates the poisson equation, thus the mass, and the
%magnification of the difflaction field created by the alphas found
%beforehand.
%load alpha_plus
%load transformed
%load transformed0201;
%load transformed10000_NEWFILE;
%load transformed2800011495FINAL
%load transformed35000_11495_10extras
load transformed10000_11495_3extras55
%load bestmodel_aCovtry
%load transformedspeedup&free10000_1902
%load transformedspeedup4free20000TRY4
%load transformedspeedup3free20000fixed1
%load transformed10000_4SS1good
%load transformedRegProc10000_3
% rat=900/alpha_x_DM(x_size-i_x-5,round((x_size-i_x)/2));
% rat_gal=10000/alpha_x(x_size-i_x-5,round((x_size-i_x)/2));
% 
% k_new=0.5
% factor=1.8*k_new;
% gamma=0.013
% phi=1.835
% %k_plus=1.18
% k_gal=0.008
% for m=2:(x_size-i_x)
%    % m
%  for  n=2:(y_size-i_y)
%      
%    alpha_x_add(m,n)=gamma*cos(2*phi)*m+gamma*sin(2*phi)*n;
%    alpha_y_add(m,n)=gamma*sin(2*phi)*m-gamma*cos(2*phi)*n;
%      
%  end
%end
% q=1
%  alpha_x_plus(1:(x_size-i_x),1:(x_size-i_x))=0;
%  alpha_y_plus(1:(y_size-i_y),1:(y_size-i_y))=0;
% for m=2:(x_size-i_x)
%    % m
%  for  n=2:(y_size-i_y)
%    
%             dif_x=-(1250-(m+i_x-1));
%             dif_y=-(1250-(n+i_y-1));
%             shoresh=sqrt(dif_x^2+dif_y^2);
%             
%             if (shoresh > 0)
%                 alpha_x_plus(m,n)=alpha_x_plus(m,n)+factor*868*(dif_x/shoresh^(1-(-q+1)));
%                 alpha_y_plus(m,n)=alpha_y_plus(m,n)+factor*868*(dif_y/shoresh^(1-(-q+1)));
%                 end
%             if (shoresh <= 0)
%                 alpha_x_plus(m,n)=alpha_x_plus(m,n)+0.0;
%                 alpha_y_plus(m,n)=alpha_y_plus(m,n)+0.0;
%                  end  
%   
%  end
% end
% alpha_x_ALL=(k_gal*alpha_x(1:2500,1:2500)*rat_gal+(1-k_gal)*alpha_x_DM(1:2500,1:2500)*rat)*factor+alpha_x_add(1:2500,1:2500);%+k_gal*k_new*k_plus*alpha_x_plus(1:2500,1:2500);
% alpha_y_ALL=(k_gal*alpha_y(1:2500,1:2500)*rat_gal+(1-k_gal)*alpha_y_DM(1:2500,1:2500)*rat)*factor+alpha_y_add(1:2500,1:2500);%+k_gal*k_new*k_plus*alpha_y_plus(1:2500,1:2500);
%this factor should be itterated. Always with plus sign: 
alpha_x_ALL=alpha_x_ALL*1.47;
alpha_y_ALL=alpha_y_ALL*1.47;
%calculate magnification and mass:
for m=2:2998
   % m
 for  n=2:2998
     %
                da_x=(alpha_x_ALL(m+1,n)-alpha_x_ALL(m-1,n))/2.0;
                da_y=(alpha_y_ALL(m,n+1)-alpha_y_ALL(m,n-1))/2.0;
                da_x_dy=(alpha_x_ALL(m,n+1)-alpha_x_ALL(m,n-1))/2.0;
                da_y_dx=(alpha_y_ALL(m+1,n)-alpha_y_ALL(m-1,n))/2.0;
                poisson_ALL(m-1,n-1)=da_x+da_y;
                magnification_ALL(m-1,n-1)=abs(1.0/(1.0-poisson_ALL(m-1,n-1)+da_x*da_y-da_x_dy*da_y_dx));
                %magnification_with_sign(m,n)=(1.0/(1.0-poisson(m,n)+da_x*da_y-da_x_dy*da_y_dx));
                 
                end
end

%changing axis just for show
%magnification_im_ALL=magnification_ALL';                
imagesc(log(magnification_ALL'))
%colormap(gray);
set(gca,'YDir','normal');
%save ('magnif_ALL_z1p89_v3FINAL.mat','magnification_ALL','alpha_x_ALL','alpha_y_ALL')
display('calculated magnification of Both. check the image.')
display('now use cluster_lense_source_to_image.m and cluster_image_to_source.m')
display('and  "cluster_image_to_source_itr_new.m" to verify images and scaling')
