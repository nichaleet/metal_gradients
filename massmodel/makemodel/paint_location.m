function paint_location
load transformedTry1BCG1RedApr16Gaussian

f=load('imagesA611.txt');
g=load('lensed_locations.txt');
st=new_im;%imread('stiff.tif');


for i=1:length(f(:,1))
    lx=(f(i,1)-i_x);
    ly=(f(i,2)-i_y);
   st(lx-10:lx+10, ly-10:ly+10,1)=0;
   st(lx-10:lx+10, ly-10:ly+10,2)=255;
   st(lx-10:lx+10, ly-10:ly+10,3)=0; 
       lx=g(i,1)-i_x;
    ly=g(i,2)-i_y;
   st(lx-10:lx+10, ly-10:ly+10,1)=255;
   st(lx-10:lx+10, ly-10:ly+10,2)=0;
   st(lx-10:lx+10, ly-10:ly+10,3)=0; 
    
end

imshow(st)