pro metallicity_map,type=type,name=name,path=path, Ha_file=Ha_file, NII_file=NII_file,Hb_file=Hb_file,OIII_file=OIII_file,Haerr_file=Haerr_file,OIIIerr_file=OIIIerr_file,obj
common BPTshare, goodagn
;From Pettini and Pagel, 2004
;INPUT
;type = 'O3N2' if want to calculate metallicity from both methods. If not, type anything
;
;note that Kn2 image is always bigger than Hbb
; E.g.
;metallicity_map,type='O3N2',name='CSWA128',path='/scr2/nichal/workspace/output/5sigmalim/1smoothparam/',Ha_file='s130911_a021_mosaic_Kn2_100_acube_obj.fits', NII_file='s130911_a021_mosaic_Kn2_100_secondline_acube_obj.fits',Hb_file='s130911_a023_mosaic_Hbb_100_secondline_acube_obj.fits',OIII_file='s130911_a023_mosaic_Hbb_100_acube_obj.fits',Haerr_file='s130911_a021_mosaic_Kn2_100_aerrorcube_obj.fits',OIIIerr_file='s130911_a023_mosaic_Hbb_100_aerrorcube_obj.fits',1


Ha_file   = path+Ha_file
NII_file  = path+NII_file
Hb_file   = path+Hb_file
OIII_file = path+OIII_file
Haerr_file   = path+Haerr_file
OIIIerr_file = path+OIIIerr_file


;============================================
;1) Metallicity from Ha and NII

hacube     = readfits(Ha_file,Ha_hdr)
NII_cube   = readfits(NII_file,NII,hdr)
haerr_cube = readfits(Haerr_file)

halpha     = hacube[ *,*,1]
halpha_err = haerr_cube[*,*,1]
NII        = NII_cube[*,*,0]
NII_err    = NII_cube[*,*,1]

N2      = ALOG10(NII/halpha)
N2_err  = 0.43429*sqrt((NII_err/NII)^2+(halpha_err/halpha)^2)
N2metal = 8.90+0.57*N2   ;metallicity unit = 12+log(O/H)
N2metal_err = 0.57*N2_err
N2fits  = [[[N2metal]],[[N2metal_err]]] 


if type ne 'O3N2' then begin
   sz = size(hacube)
   sizex = sz(1)
   sizey = sz(2)
   data_quality, path=path,fits_file=[ha_file,nii_file,'no','no'],[1,2,0,0],[0,0,0,0],[0,0,0,0],sizex,sizey,nameplot=[name+' Halpha',name+' NII',name+' OIII',name+' Hbeta'],namegal=name
endif


PA_spec    = sxpar(ha_hdr,'PA_SPEC')
scale      = sxpar(ha_hdr,'scale') ; in degrees
scale      = scale*3600. ; in arcsecs

RA_Ha_ref  = sxpar(ha_hdr,'CRVAL2') ;R.A. at reference pixel in mosaic frame in degrees
DEC_Ha_ref = sxpar(ha_hdr,'CRVAL3') ;DEC at reference pixel in mosaic frame in degress
x_ha_ref   = sxpar(ha_hdr,'CRPIX2') ;Reference pixel location     
y_ha_ref   = sxpar(ha_hdr,'CRPIX3') ;Reference pixel location     

;=============================================
;2) Metallicity from Ha, NII, Hb, and OIII
if type eq "O3N2" then begin
hbcube = readfits(Hb_file,Hb_hdr)
OIIIcube = readfits(OIII_file,OIII_hdr)
OIIIerr_cube = readfits(OIIIerr_file)

OIII = OIIIcube[*,*,1]
OIII_err = OIIIerr_cube[*,*,1]
hbeta = hbcube[*,*,0]
hbeta_err = hbcube[*,*,1]


;--------------------------------------------
;Align pixels from 2 sets

if sxpar(ha_hdr,'PA_SPEC') ne sxpar(hb_hdr,'PA_SPEC') or sxpar(ha_hdr,'PA_IMAG') ne sxpar(hb_hdr,'PA_IMAG') then begin
   print, 'Error: POSITION ANGLES OF 2 IMAGES DO NOT ALIGN'
   goto, skip
endif

if sxpar(ha_hdr,'scale') ne sxpar(hb_hdr,'scale') then begin
   print, 'Error: 2 Images do not have the same scale.'
   goto, skip
endif


RA_Hb_ref  = sxpar(hb_hdr,'CRVAL2') ;R.A. at reference pixel in mosaic frame in degrees
DEC_Hb_ref = sxpar(hb_hdr,'CRVAL3') ;DEC at reference pixel in mosaic frame in degrees
x_hb_ref   = sxpar(hb_hdr,'CRPIX2') ;Reference pixel location     
y_hb_ref   = sxpar(hb_hdr,'CRPIX3') ;Reference pixel location     

;;;;;;;;;;;;;;;;;;;;
;If the input image is obj (cut off version from mosaic then we have to correct its reference y position
if obj eq 1 then begin
   pix_y1_ha = sxpar(ha_hdr,'PIX_Y1')
   pix_y1_hb = sxpar(hb_hdr,'PIX_Y1')
   y_ha_ref = y_ha_ref-pix_y1_ha
   y_hb_ref = y_hb_ref-pix_y1_hb
endif
;;;;;;;;;;;;;;;;;;;
;Find pixel position of RA_Hb_ref, DEC_Hb_ref in Halpha image

GCIRC,1,RA_Ha_ref/15.,Dec_Ha_ref,RA_Hb_ref/15.,Dec_Hb_ref,d_hbref_haref ; d is in arcsec
POSANG,1,RA_Ha_ref/15.,Dec_Ha_ref,RA_Hb_ref/15.,Dec_Hb_ref,theta_prime

theta = theta_prime-PA_spec  ; in degrees
x_Hbref_inHa  = d_hbref_haref*cos(theta*!dpi/180.)/scale+x_Ha_ref
y_Hbref_inHa  = d_hbref_haref*sin(theta*!dpi/180.)/scale+y_ha_ref

print,'(x,y) in Hbeta:', x_hb_ref,y_hb_ref
print,'matched(x,y) in Halpha', x_Hbref_inHa,y_Hbref_inHa

x_hbref_inHa = round(x_hbref_inHa)
y_hbref_inHa = round(y_hbref_inHa)

; Make an array of new intersected image
size_ha = size(halpha)
size_hb = size(hbeta)

x_lat_ha = size_ha(1)-x_hbref_inHa
x_lat_hb = size_hb(1)-x_hb_ref
size_x = min([x_hbref_inHa,x_hb_ref])+min([x_lat_ha,x_lat_hb])

y_lat_ha = size_ha(2)-y_hbref_inHa
y_lat_hb = size_hb(2)-y_hb_ref
size_y = min([y_hbref_inHa,y_hb_ref])+min([y_lat_ha,y_lat_hb])

;x_arr = rebin(findgen(size_x),size_x,size_y)
;y_arr = transpose(rebin(findgen(size_y),size_y,size_x))

x_new_ref =  min([x_hbref_inHa,x_hb_ref])
y_new_ref =  min([y_hbref_inHa,y_hb_ref])

x_shift_ha = x_hbref_inHa-x_new_ref
y_shift_ha = y_hbref_inHa-y_new_ref
x_shift_hb = x_hb_ref-x_new_ref
y_shift_hb = y_hb_ref-y_new_ref

halpha_new     = findgen(size_x,size_y)
halpha_err_new = halpha_new
hbeta_new      = findgen(size_x,size_y)
hbeta_err_new  = hbeta_new
OIII_new       = findgen(size_x,size_y)
OIII_err_new   = OIII_new
NII_new        = findgen(size_x,size_y)
NII_err_new    = NII_new
N2fits_new     = findgen(size_x,size_y,2)
N2metal_new    = findgen(size_x,size_y)
N2metal_err_new= N2metal_new

for i=0,size_x-1 do begin for j=0,size_y-1 do begin

halpha_new(i,j) = halpha(i+x_shift_ha,j+y_shift_ha)
halpha_err_new(i,j) = halpha_err(i+x_shift_ha,j+y_shift_ha)

NII_new(i,j)    = NII(i+x_shift_ha,j+y_shift_ha)
NII_err_new(i,j)    = NII_err(i+x_shift_ha,j+y_shift_ha)

N2fits_new(i,j,0) = N2fits(i+x_shift_ha,j+y_shift_ha,0)
N2fits_new(i,j,1) = N2fits(i+x_shift_ha,j+y_shift_ha,1)
N2metal_new(i,j)      = N2metal(i+x_shift_ha,j+y_shift_ha)
N2metal_err_new(i,j)  = N2metal_err(i+x_shift_ha,j+y_shift_ha)

hbeta_new(i,j)  = hbeta(i+x_shift_hb,j+y_shift_hb)
hbeta_err_new(i,j)  = hbeta_err(i+x_shift_hb,j+y_shift_hb)

OIII_new(i,j)   = OIII(i+x_shift_hb,j+y_shift_hb)
OIII_err_new(i,j)   = OIII_err(i+x_shift_hb,j+y_shift_hb)

endfor
endfor

NII_err_new(where(NII_New eq 0.)) = 1./0.
NII_new(where(NII_New eq 0.)) = 1./0.
Hbeta_err_new(where(hbeta_new eq 0.)) = 1./0.
Hbeta_new(where(Hbeta_new eq 0.)) = 1./0.


;--------------------------------------------
;CALCULATE METALICITY
O3N2 = alog10((OIII_new/hbeta_new)/(NII_new/halpha_new))
O3N2metal = 8.73-0.32*O3N2

O3N2_err = 0.43429*sqrt((OIII_err_new/OIII_new)^2+(hbeta_err_new/hbeta_new)^2+(NII_err_new/NII_new)^2+(halpha_err_new/halpha_new)^2)
O3N2metal_err = 0.32*O3N2_err
 
;; difference in metallicity using 2 methods
Diff_metal = O3N2metal-N2metal_new
Diff_metal_err = sqrt(O3N2metal_err^2+N2metal_err_new^2)

;;Halpha/Hbeta map
Ha_over_Hb = Halpha_new/Hbeta_new
Ha_over_Hb_err = sqrt((Halpha_err_new/Halpha)^2+(Hbeta_err_new/Hbeta)^2)


; MAKING BPT DIAGRAMS

N2divHa =  NII_new/halpha_new
O3divHb =  OIII_new/hbeta_new
N2divHa_err = sqrt((NII_err_new/NII_new)^2+(halpha_err_new/halpha_new)^2)*N2divHa
O3divHb_err = sqrt((OIII_err_new/OIII_new)^2+(hbeta_err_new/hbeta_new)^2)*O3divHb
pix_indices = indgen((size(nii_new))[1]*(size(nii_new))[2])

nan_ind1 = where(~finite(N2divHa))
nan_ind2 = where(~finite(O3divHb))
nan_ind = setunion(nan_ind1,nan_ind2)
if n_elements(nan_ind) ne n_elements(n2divha) then remove, nan_ind,N2divHa,O3divHb,N2divHa_err,O3divHb_err,pix_indices else begin 
   n2divha = []
   o3divhb = []
endelse

; NOW ALL n2divha,o3divhb are not NAN but have smaller sizes. Their original indices in the metallicity map are stored in pix_indices

if n_elements(N2divHa) le 1. then begin
   print, '!!!!'+name+' does not have enough good data to plot BPT diagram'
   agnflag = 0
   map_agn = fltarr((size(nii_new))[1],(size(nii_new))[2])
   map_agn(*,*) = 1./0.
   goto, nobptdata
endif

  ; Find points that are above Kewley AGN line
AGN = 0.61/(alog10(n2divha)-0.47) + 1.19
agnpointsinn2divha = (where(alog10(O3divhb) ge agn)) ;indices in n2divha

if agnpointsinn2divha ne [-1] then begin   
   AGNpoints = pix_indices(where(alog10(O3divhb) ge agn)) ;indices in metal map
   map_agn = fltarr((size(nii_new))[1],(size(nii_new))[2])
   map_agn(AGNpoints) = 1.
   agnflag= 1
endif else begin 
   print, 'NO AGN. Yay!'
   map_agn = fltarr((size(nii_new))[1],(size(nii_new))[2])
   agnflag = 0
endelse

; Analyze indicated AGN points and make intensity vs wavelengths and fitted gaussian profiles
if agnflag eq 1 then begin
   agnpointsxy = array_indices(map_agn,agnpoints)  ;indices in metal map (x,y)
   agn_fitplots,  path=path,fits_file=[ha_file,nii_file,oiii_file,hb_file],[1,2,1,2],[x_shift_ha,x_shift_ha,x_shift_hb,x_shift_hb],[y_shift_ha,y_shift_ha,y_shift_hb,y_shift_hb],size_x,size_y,nameplot=[name+' Halpha',name+' NII',name+' OIII',name+' Hbeta'],namegal=name,agnpointsxy,N2divHa(agnpointsinn2divha),O3divHb(agnpointsinn2divha),N2divHa_err(agnpointsinn2divha),O3divHb_err(agnpointsinn2divha),agnpointsinn2divha

goodagninn2divha = agnpointsinn2divha(goodagn)
endif
;stop
; PLOT BPT DIAGRAM

set_plot, 'ps'
device, file='output/metallicity/BPT_'+name+'.eps', encapsulated =1,/color
device, xsize=15, ysize=20
!p.multi = [0,1,2,0,0]
!p.font = 0
!x.margin = [12,3]
!p.thick = 4
!x.thick = 2
!y.thick = 2

rangex = [min(n2divha(where(n2divha gt 0)))*0.9,max(n2divha(where(n2divha gt 0)))*1.2]
rangey = [min(o3divhb(where(o3divhb gt 0)))*0.9,max(o3divhb(where(o3divhb gt 0)))*1.2]

plot,N2divHa,O3divHb,psym=1,xtitle='[NII/Ha]',ytitle='[OIII/Hb]',title=name,/xlog,/ylog,xrange=rangex,yrange=rangey

if agnflag eq 1 then oploterror,N2divHa(goodagninn2divha),O3divHb(goodagninn2divha),N2divHa_err(goodagninn2divha),o3divHb_err(goodagninn2divha),psym=4,thick=5

  ; Show the division between AGN and starbursts.
   ; Kauffmann et al (2003) line:
n2ha = findgen(150)/100.
log_o3hb_kauff = 0.61/(alog10(n2ha)-0.05) + 1.3
oplot, n2ha, 10.^log_o3hb_kauff, linestyle=2,min_value=!y.crange[0]
   ; Kewley et al (2001) line:
n2ha = findgen(200)/100.
log_o3hb = 0.61/(alog10(n2ha)-0.47) + 1.19
oplot, n2ha, 10.^log_o3hb, linestyle=1
legend,['Kauffmann','Kewley'],linestyle=[2,1]

ploterror,N2divHa,O3divHb,N2divHa_err,O3divHb_err,/xlog,/ylog,psym=1,xrange=rangex,yrange=rangey,xtitle='[NII/Ha]',ytitle='[OIII/Hb]',title=name ;type=3 is log log scale

if agnflag eq 1 then begin 
   oplot,n2divha(goodagninn2divha),o3divhb(goodagninn2divha),psym=4,thick=5 ;make significant AGN points diamond in BPT diagram
   print, 'Coordinates of diamond points in BPT diagram (x,y):'
   print, agnpointsxy(*,goodagn)
;stop
   
endif
oplot, n2ha, 10^log_o3hb_kauff,linestyle=2,min_value=rangey[0]
oplot, n2ha, 10^log_o3hb,linestyle=1
if agnflag eq 1 then xyouts,2,5,'Total number of AGN pixels with significances of  NII and Hb fits >0 = '+string(n_elements(goodagn),format='(I2)'),/device
device,/close

set_plot,'x'
plot,N2divHa,O3divHb,psym=1,background=getcolor('white'),color=0,xrange=[0,1.5],yrange=[0,20],xtitle='NII/Ha',ytitle='OIII/Hb',title=name
;ploterror,N2divHa,O3divHb,N2divHa_err,O3divHb_err,psym=1,background=getcolor('white'),color=0,xtitle='NII/Ha',ytitle='OIII/Hb',title=name


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;analyze data quality: histograms of fit significance/ sum of data and fit for each line

nobptdata:
 
data_quality, path=path,fits_file=[ha_file,nii_file,oiii_file,hb_file],[1,2,1,2],[x_shift_ha,x_shift_ha,x_shift_hb,x_shift_hb],[y_shift_ha,y_shift_ha,y_shift_hb,y_shift_hb],size_x,size_y,nameplot=[name+' Halpha',name+' NII',name+' OIII',name+' Hbeta'],namegal=name



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Creating fits array
O3N2_fits = [[[N2metal_new]],[[O3N2metal]],[[Diff_metal]],[[Ha_over_hb]],[[map_agn]]]
O3N2_fits_err = [[[N2metal_err_new]],[[O3N2metal_err]],[[Diff_metal_err]],[[Ha_over_Hb_err]]]

;MAKE HEADER
mkhdr,o3n2_header,O3N2_Fits
mkhdr,O3N2_err_header,O3N2_fits_err
hdr_para_names = ['frame1','frame2','frame3','frame4','frame5','SCALE','REF_x','REF_y','REF_RA','REF_DEC','PA_SPEC']
hdr_para_vals = [1/0,1/0,1/0,1/0,1/0,scale,x_new_ref,y_new_ref,RA_Hb_ref,DEC_Hb_ref,pa_spec]
comment1 = ['N2metal','O3N2metal','O3N2metal-N2metal','Ha/Hb','Detected AGN points','arcsecs per pixel','position of reference pixel in this image','position of reference pixel in this image','RA of the above reference pixel','DEC of the above reference pixel','position angle of spectograph and this image']
comment2 = ['N2metal uncertainty','O3N2metal uncertainty','O3N2metal-N2metal uncertainty','Ha/Hb uncertainty','nothing','arcsecs per pixel','position of reference pixel in this image','position of reference pixel in this image','RA of the above reference pixel','DEC of the above reference pixel','position angle of spectograph and this image']
for ii=0,n_elements(hdr_para_names)-1 do begin 
   sxaddpar,o3n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment1(ii)
   sxaddpar,o3n2_err_header,hdr_para_names(ii),hdr_para_vals(ii),comment2(ii)
endfor


writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicity.fits',O3N2_fits,o3n2_header
writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicity_error.fits',O3N2_fits_err,o3n2_err_header


;;;;;;;;;;;;;;;;;;
;Make Metallicity as a function of radius
if name eq 'CSWA128' then xy_Hacont_peak = [22,35] else $
if name eq 'CSWA165' then xy_Hacont_peak = [25,37] else goto, skip

ypeak_pos = xy_Hacont_peak(1)-pix_y1_ha-y_shift_ha
xpeak_pos = xy_hacont_peak(0)-x_shift_ha
print, 'xy position of Ha cont peak in O3N2_metallicity.fits', xpeak_pos, ypeak_pos
x_array = (findgen(size_x)-xpeak_pos)/10. ; arcsecs from Ha cont peak position

metallicity_rad = findgen(size_x)
metallicity_rad_err = findgen(size_x)
n_average = findgen(size_x)

for ii=0,size_x-1 do begin
   metallicity_area = findgen(5)
   metallicity_area_err = findgen(5)
   for jj=0,4 do begin
      metallicity_area(jj) = N2metal_new(ii,ypeak_pos+jj-2)
      metallicity_area_err(jj) = N2metal_err_new(ii,ypeak_pos+jj-2)
   endfor
   bad_pix = where(~finite(metallicity_area))
   if n_elements(bad_pix) eq 5 then begin
      metallicity_rad(ii) = 0.
      metallicity_rad_err(ii) = 0.
      n_average(ii) = 0
   endif else begin
      if bad_pix(0) ne -1 then remove, bad_pix,metallicity_area,metallicity_area_err
      ;metallicity_rad(ii)  = mean(metallicity_area)
      meanerr,metallicity_area,metallicity_area_err,mean_metal,sigma_mean_metal,sigma_Data
      metallicity_rad(ii) = mean_metal
      metallicity_rad_err(ii) = sigma_data
      n_average(ii) = n_elements(metallicity_area)
   endelse
endfor

o3metallicity_rad = findgen(size_x)
o3metallicity_rad_err = findgen(size_x)
o3n_average = findgen(size_x)

for ii=0,size_x-1 do begin
   metallicity_area = findgen(5)
   metallicity_area_err = findgen(5)
   for jj=0,4 do begin
      metallicity_area(jj) = O3N2metal(ii,ypeak_pos+jj-2)
      metallicity_area_err(jj) = O3N2metal_err(ii,ypeak_pos+jj-2)
   endfor
   bad_pix = where(~finite(metallicity_area))
   if n_elements(bad_pix) eq 5 then begin
      o3metallicity_rad(ii) = 0.
      o3metallicity_rad_err(ii) = 0.
      o3n_average(ii) = 0
   endif else begin
      if bad_pix(0) ne -1 then remove, bad_pix,metallicity_area,metallicity_area_err
      ;metallicity_rad(ii)  = mean(metallicity_area)
      meanerr,metallicity_area,metallicity_area_err,mean_metal,sigma_mean_metal,sigma_Data
      o3metallicity_rad(ii) = mean_metal
      o3metallicity_rad_err(ii) = sigma_data
      o3n_average(ii) = n_elements(metallicity_area)
   endelse
endfor


set_plot, 'ps'
device, file='output/metallicity/Metallicity_gradient_'+name+'.eps', encapsulated =1, /color
device, xsize=30, ysize=20
!p.multi = [0,2,2,0,0]
plot, x_array, metallicity_rad,psym=1,xtitle='distance from Ha Cont peak(arcsecs)', ytitle = 'N2 metallicity', title = name,yrange=[7,9],xrange=[-2,4]
ploterror,x_array,metallicity_rad,metallicity_rad_err,psym=1,xtitle='distance from Ha Cont peak(arcsecs)', ytitle = 'N2 metallicity', title = name,yrange=[7,9],xrange=[-2,4]
oplot,x_array,n_average/10.+7.

plot, x_array, o3metallicity_rad,psym=1,xtitle='distance from Ha Cont peak(arcsecs)', ytitle = 'O3N2 metallicity',yrange=[7,9],xrange=[-2,4]
ploterror,x_array,o3metallicity_rad,o3metallicity_rad_err,psym=1,xtitle='distance from Ha Cont peak(arcsecs)', ytitle = 'O3N2 metallicity',yrange=[7,9],xrange=[-2,4]
oplot,x_array,o3n_average/10.+7.

device,/close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endif 

skip:


if type ne "O3N2" then begin
   hdr_para_vals = [scale,x_ha_ref,y_ha_ref,RA_HA_ref,DEC_Ha_ref,pa_spec] 
   mkhdr,n2_header,n2fits
   hdr_para_names = ['SCALE','REF_x','REF_y','REF_RA','REF_DEC','PA_SPEC']
   comment = ['arcsecs per pixel','position of reference pixel in this image','position of reference pixel in this image','RA of the above reference pixel','DEC of the above reference pixel','position angle of spectograph and this image']

   for ii=0,n_elements(hdr_para_names)-1 do sxaddpar,n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment(ii)

   writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_N2_metallicity.fits',N2fits,n2_header
endif

end



