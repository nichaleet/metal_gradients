function mark_polygons_of_images(allimagesfile,imtomark)


load transformedTry1BCG1RedApr15
f=load(allimagesfile);
%g=imread(imtomark);

no_of_sys=f(1,5);
line_no=1;
for i=1:no_of_sys
    no_of_im_insys=f(line_no,4);
    BWnew(1:leng,1:leng,:)=0;
    for j=1:no_of_im_insys
        line_no
xc=f(line_no,1)-x_start;
yc=f(line_no,2)-y_start;
sizec=200;
display('mark the image');
try_im=new_im(xc-round(sizec/2):xc+round(sizec/2),yc-round(sizec/2):yc+round(sizec/2),:);
BW=roipolyold(try_im); 
display('OK')
BWnew(xc-round(sizec/2):xc+round(sizec/2),yc-round(sizec/2):yc+round(sizec/2),:)=BW(:,:,:);
line_no=line_no+1;
    end
    %save matrix BWnew for each system and change index
    save('BWnew1.mat','BWnew') 
end