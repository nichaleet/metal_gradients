function transform_image(image_name);
%load bestmodel_aCovtry
load bestmodel_aCovtry_newallsysV1short.mat
%load model_example_gaussian0411
clear alpha_x_DM alpha_x_add alpha_y_DM alpha_y_add X Y XI YI da_x_dx da_y_dy da_x_dy da_y_dx dif_x dif_y poisson
clear x_factor y_factor x_factor_shift y_factor_shift
img=imread(image_name);

temp_im=img;

for m=1:x_size-i_x
 for n=1:y_size-i_y
   new_im(m,n,:)=temp_im((leng+2)-n,m,:);
     
 end
end
save transformed0201;
display('now you can use "lensing_tool.m"')
