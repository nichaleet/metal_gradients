pro mosaic_mac1133
;Below is for pipeline mosaic files
File1 = "/scr2/nichal/workspace/reduced_data/Feb2014/scaledskysubtracted/s140223_a016_mosaic_Hn4_100.fits"
file2 = "/scr2/nichal/workspace/reduced_data/Mar2013/scaledskysubtracted/s130303_a016_mosaic_Hn4_100.fits"


cube1 = readfits(file1,head1)
cube2 = readfits(file2,head2)

;get header parameters

RA_1     = sxpar(head1,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC_1    = sxpar(head1,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
PA_Spec_1 = sxpar(head1,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG_1 = sxpar(head1,'PA_IMAG') ; position angle of imager on sky
NAXIS_1   = sxpar(head1,'NAXIS*')

CRVAL1_1 = sxpar(head1, 'CRVAL1') ;[nm] Wavelength at reference pixel 
CRVAL2_1 = sxpar(head1, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3_1 = sxpar(head1, 'CRVAL3') ;[deg] DEC at reference pixel       

CRPIX1_1 = sxpar(head1, 'CRPIX1') ; Reference pixel location        
CRPIX2_1 = sxpar(head1, 'CRPIX2') ; Reference pixel location        
CRPIX3_1 = sxpar(head1, 'CRPIX3')  ;Reference pixel location       

CDELT1_1 = sxpar(head1, 'CDELT1') ; Wavelength scale in nm/channel
CDELT2_1 = sxpar(head1, 'CDELT2') ; Pixel scale in degree/pixel
CDELT3_1 = sxpar(head1, 'CDELT3') ; Pixel scale in degree/pixel
PC2_2_1    =  sxpar(head1, 'PC2_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC2_3_1    =  sxpar(head1, 'PC2_3')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_2_1    =  sxpar(head1, 'PC3_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_3_1    =  sxpar(head1, 'PC3_3')  ;/RA, Dec axes rotated by 205.000 degr.    
angle1 = (205.-90.)/180.*!dpi       


RA_2     = sxpar(head2,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC_2    = sxpar(head2,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
PA_Spec_2 = sxpar(head2,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG_2 = sxpar(head2,'PA_IMAG') ; position angle of imager on sky
NAXIS_2   = sxpar(head2,'NAXIS*')

CRVAL1_2 = sxpar(head2, 'CRVAL1') ;[nm] Wavelength at reference pixel 
CRVAL2_2 = sxpar(head2, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3_2 = sxpar(head2, 'CRVAL3') ;[deg] DEC at reference pixel       

CRPIX1_2 = sxpar(head2, 'CRPIX1') ; Reference pixel location        
CRPIX2_2 = sxpar(head2, 'CRPIX2') ; Reference pixel location        
CRPIX3_2 = sxpar(head2, 'CRPIX3')  ;Reference pixel location       

CDELT1_2 = sxpar(head2, 'CDELT1') ; Wavelength scale in nm/channel
CDELT2_2 = sxpar(head2, 'CDELT2') ; Pixel scale in degree/pixel
CDELT3_2 = sxpar(head2, 'CDELT3') ; Pixel scale in degree/pixel
PC2_2_2   =  sxpar(head2, 'PC2_2')  ;/RA, Dec axes rotated by  degr.           
PC2_3_2   =  sxpar(head2, 'PC2_3')  ;/RA, Dec axes rotated by 25.000 degr.           
PC3_2_2   =  sxpar(head2, 'PC3_2')  ;/RA, Dec axes rotated by 25.000 degr.           
PC3_3_2   =  sxpar(head2, 'PC3_3')  ;/RA, Dec axes rotated by 25.000 degr.      
angle2 = (25.-90.)/180.*!dpi     


;From p.43 in research log, found that the difference in PA's is 180 degree -> The image of Feb is the upside down of Mar. So we just need to find a matching pixel then align them (reversely: +x is -x in another picture, +y is -y in another picture)
scalex_1 = 0.1
scaley_1 = 0.1
pixmatch,pa_spec_1,scalex_1,scaley_1,CRVAL2_1,CRVAL2_2,CRVAL3_1,CRVAL3_2,CRPIX2_1,CRPIX3_1,CRPIX2_2,CRPIX3_2,x_of_ref2_in1,y_of_ref2_in1
print,'Feb file has RA,DEC: ',CRVAL2_1,CRVAL3_1
print,'Mar file has RA,DEC: ',CRVAL2_2,CRVAL3_2

;To plot contour map
ha_sum1 = total(cube1[89:95,*,*],1)
ha_sum2 = total(cube2[89:95,*,*],1)
mkhdr,ha_sum1hdr,ha_sum1
mkhdr,ha_sum2hdr,ha_sum2
hdr_paranames=['CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','PC1_1','PC1_2','PC2_1','PC2_2']
;hdr_val1 = [CDELT2_1,CDELT3_1,CRPIX2_1,CRPIX3_1,CRVAL2_1,CRVAL3_1,PC2_2_1,PC2_3_1,PC3_2_1,PC3_3_1]
hdr_val1 = [CDELT2_1,CDELT3_1,CRPIX2_1,CRPIX3_1,CRVAL2_1,CRVAL3_1,cos(angle1),sin(angle1),-1.*sin(angle1),cos(angle1)]

;hdr_val2 = [CDELT2_2,CDELT3_2,43,64,CRVAL2_1,CRVAL3_1,PC2_2_2,PC2_3_2,PC3_2_2,PC3_3_2]
hdr_val2 = [CDELT2_2,CDELT3_2,43,64,CRVAL2_1,CRVAL3_1,cos(angle2),sin(angle2),-1.*sin(angle2),cos(angle2)]
hdr_paraname_str = ['CTYPE1','CTYPE2','CUNIT1','CUNIT2']
hdr_paraname_str_val = ['RA---TAN','DEC--TAN','deg','deg']

for ii=0,n_elements(hdr_paranames)-1 do begin 
   sxaddpar,ha_sum1hdr,hdr_paranames(ii),hdr_val1(ii)
   sxaddpar,ha_sum2hdr,hdr_paranames(ii),hdr_val2(ii)
endfor
for ii=0,n_elements(hdr_paraname_Str)-1 do begin 
   sxaddpar,ha_sum1hdr,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
   sxaddpar,ha_sum2hdr,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
endfor
writefits,'FEbmap.fits',ha_sum1,ha_sum1hdr

writefits, 'Marchmap.fits',ha_sum2,ha_sum2hdr

;Adding the 2 cubes together into a new mosaic-ed cube
sz = size(cube1)
newcube = fltarr(sz[1],58,61)
for i=0,57 do for j=0,60 do begin
newcube(*,i,j)=cube1(*,8+i,16+j)+cube2(*,67-i,76-j)
endfor
writefits,'macs1133_Ha_Hn4_handmosaic_scaledsky_and_simplyaddition_4hr.fits',newcube,head1

;plot contour map
ha_sum = total(newcube[90:94,*,*],1)
imdisp,ha_sum,out_pos=pos1
;contour, ha_sum, pos=pos1, /xs, /ys, /noerase, color = 255
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.02,0.03,0.04,0.05],color=255
writefits,'newcubeHa.fits',ha_sum

stop

end
