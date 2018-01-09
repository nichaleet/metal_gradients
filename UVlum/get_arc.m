
function get_arc(image_name,name,i_x,f_x,i_y,f_y);
% this program will show the image, let you mark the image you are
% interested in 
img=fitsread(image_name,'image',1);
img = flipud(img);
[y_size,x_size]=size(img);
x_leng = f_x-i_x+1;
y_leng = f_y-i_y+1;
final_roi = zeros(x_size,y_size);
for m=1:y_leng
   for n=1:x_leng
     im(m,n)=img(y_size-m-i_y,n+i_x);
   end
end

imbit = (im-min(min(im)))/max(max(im));
imbit = exp(imbit); %separate the background from signal
imbit = (imbit-min(min(imbit)))/max(max(imbit))*1000.;
imbit = fix(imbit) ;% now background will be zero but the signal has been exponentialed. So the bright spot has really high value while the faint spot has a low value.

imbit(imbit<4) = NaN;
imbit = imbit+2.   ;%The lowest was 1 so if we take log, it would be zero. so had to add 1 so thatthe lowest is now 3.
imbit = log(imbit);
imbit(isnan(imbit))=0.;
imbit = (imbit-min(min(imbit)))/max(max(imbit))*255;

imshow(imbit);
file_roi = strcat(name,'_roi.mat')
if exist(file_roi,'file')
    load (file_roi)
else
    BW=roipoly(imbit);
    BW = double(BW);
end

for m=1:y_leng
    for n=1:x_leng
        final_roi(i_x+n-1,i_y+m-1)=BW(m,n);
    end
end
save(file_roi,'final_roi');
fitswrite(flipud(imrotate(final_roi,90)),strcat(name,'_roi.fits'))




