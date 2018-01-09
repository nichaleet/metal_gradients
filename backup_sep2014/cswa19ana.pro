pro cswa19ana
setplot, 14
!p.multi = [0,2,2]
z=2.30
;1) do H alpha image
;loadct, 22,ncolors=100

image = readfits("/scr2/nichal/workspace/reduced_data/mosaic/cswa19_Ha_Kn1_mosaic_scaledsky_4hr.fits",head)


;get header parameters

RA     = sxpar(head,'RA')  ;RA at spatial [0,0] in mosaic ;format='(f14.10)'
DEC    = sxpar(head,'DEC');DEC at spatial [0,0] in mosaic ;format='(f14.10)'
PA_Spec = sxpar(head,'PA_SPEC') ;position angle of spectrograph on sky
PA_IMAG = sxpar(head,'PA_IMAG') ; position angle of imager on sky
NAXIS   = sxpar(head,'NAXIS*')

CRVAL1 = sxpar(head, 'CRVAL1') ;[nm] Wavelength at reference pixel 
CRVAL2 = sxpar(head, 'CRVAL2') ;[deg] R.A. at reference pixel       
CRVAL3 = sxpar(head, 'CRVAL3') ;[deg] DEC at reference pixel       

CRPIX1 = sxpar(head, 'CRPIX1') ; Reference pixel location        
CRPIX2 = sxpar(head, 'CRPIX2') ; Reference pixel location        
CRPIX3 = sxpar(head, 'CRPIX3')  ;Reference pixel location       

CDELT1 = sxpar(head, 'CDELT1') ; Wavelength scale in nm/channel
CDELT2 = sxpar(head, 'CDELT2') ; Pixel scale in degree/pixel
CDELT3 = sxpar(head, 'CDELT3') ; Pixel scale in degree/pixel
PC2_2   =  sxpar(head, 'PC2_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC2_3    =  sxpar(head, 'PC2_3')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_2    =  sxpar(head, 'PC3_2')  ;/RA, Dec axes rotated by 205.000 degr.           
PC3_3   =  sxpar(head, 'PC3_3')  ;/RA, Dec axes rotated by 205.000 degr.           

for i =0, n_elements(image[*,0,0])-1 do image[i,*,*] = filter_image(reform(image[i,*,*]),fwhm_gaussian=3)
;To plot contour. Take the Ha wavelengths
ha_sum = total(image[140:151,*,*],1)
;imdisp, ha_sum,position=[.05,.05,.5,1.], out_pos=pos1
imdisp, ha_sum, out_pos=pos1
contour, ha_sum, pos=pos1, /xs, /ys, /noerase, levels=[0.03,0.04,0.05,0.08]

;make fit file of the ha_map
ha_sum = transpose(ha_sum)
mkhdr,ha_sum_hdr,ha_sum
hdr_paranames=['CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','PC1_1','PC1_2','PC2_1','PC2_2']
hdr_val1 = [CDELT3,CDELT2,CRPIX3,CRPIX2,CRVAL2,CRVAL3,PC2_2,PC2_3,PC3_2,PC3_3]
hdr_paraname_str = ['CTYPE1','CTYPE2','CUNIT1','CUNIT2']
hdr_paraname_str_val = ['RA---TAN','DEC--TAN','deg','deg']

for ii=0,n_elements(hdr_paranames)-1 do begin 
   sxaddpar,ha_sum_hdr,hdr_paranames(ii),hdr_val1(ii)
endfor
for ii=0,n_elements(hdr_paraname_Str)-1 do begin 
   sxaddpar,ha_sum_hdr,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
endfor
writefits,'cswa19.fits',ha_sum,ha_sum_hdr


;To plot wavelengths vs intensity map. Take only galaxy region. 
ha_spectrum = fltarr(n_elements(image[*,0,0]))
for i =0, n_elements(image[*,0,0])-1 do ha_spectrum[i] = total(image[i,25:45,25:35])
wavelength = findgen(401)*0.25+1955
plot,wavelength[120:180],ha_spectrum[120:180],psym=1,/noerase

;find the peak of the ha_spectrum
wavepeak = wavelength(where(ha_spectrum eq max(ha_spectrum)))
z= wavepeak/656.4-1.
vline,(z[0]+1.)*[656.4,658.5]
xyouts,1700,1.,'MACS1133'
Print, 'redshift is ',z
print, where(ha_spectrum eq max(ha_spectrum))
stop
;wait, 30

;wait, 30


end
