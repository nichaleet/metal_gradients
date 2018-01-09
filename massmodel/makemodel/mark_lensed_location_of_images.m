function mark_lensed_location_of_images(allimagesfile,imtomark)


%load transformedTry1BCG1RedApr15
load transformed
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
line_no=line_no+1;
    end

end