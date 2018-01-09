function transform_imagefromMCspeedup(image_name);

%load BestModel7Free4RerdsFree20000TRY4
%load bestModelRegProc10000_3tryignore
%load bestMode30000_fixed1
%load bestModel28000FINAL11495
%load bestModel35000_11495_10extras
load bestModel
% [da_x_dy,da_x_dx]=gradient(alpha_x_ALL);
% [da_y_dy,da_y_dx]=gradient(alpha_y_ALL);
% poisson_ALL=da_x_dx+da_y_dy;
% magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
img=imread(image_name);

temp_im=img;

for m=1:x_size-i_x
 for n=1:y_size-i_y
   new_im(m,n,:)=temp_im((leng+2)-n,m,:);
     
 end
end
%save transformedspeedup1002;
save transformedLTMv18102013
%save transformedTry1BCG1RedApr20Gaussiancont3
display('now you can use "lensing_tool.m"')
