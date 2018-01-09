pro metallicity_map,type=type,name=name,path=path, Ha_file=Ha_file, OIII_file=OIII_file
common BPTshare, goodagn
;From Pettini and Pagel, 2004
;INPUT
;type = 'O3N2' if want to calculate metallicity from both methods. If not, type anything
;
;note that Kn2 image is always bigger than Hbb
; E.g.
;metallicity_map, type="O3N2",name='CSWA19',path='/scr2/nichal/workspace/output/', Ha_file='cswa19_Ha_Kn1_mosaic_scaledsky_4hr_acube.fits',oiii_file='cswa19_Hb_Hn1_mosaic_sky_1hr_acube.fits'

;OUTPUT:
;For O3N2 type: 
;1)fits file of  ['N2metal','O3N2metal','O3N2metal-N2metal','Ha/Hb','Detected AGN points','N2metal Error','O3N2metal error','O3N2metal-N2metal error','Ha/Hb error']
;2)Oiii_file_aligned.fits - the input oiii file but with the new size (number of pixels that each pixel align with the input Ha file)

Ha_file   = path+Ha_file
OIII_file = path+OIII_file


;============================================
;1) Metallicity from Ha and NII

hacube     = readfits(Ha_file,Ha_hdr)

halpha     = hacube[ *,*,1]
halpha_err = hacube[*,*,9]
NII        = hacube[*,*,6]
NII_err    = hacube[*,*,11]

;make a map with NII with upper limit
NII_upper = NII_err+NII
change2upper = where(NII lt 0. and NII_upper gt 0.)
NII(change2upper) = NII_upper(change2upper)
upper_tags = halpha*0.
upper_tags(where(finite(upper_tags))) = 1.
upper_tags(change2upper) = 2.
nan = where(NII_upper lt 0.)
upper_tags(nan) = 3.

N2      = ALOG10(NII/halpha)
N2_err  = 0.43429*sqrt((NII_err/NII)^2+(halpha_err/halpha)^2)
N2metal = 8.90+0.57*N2   ;metallicity unit = 12+log(O/H)
N2metal_err = 0.57*N2_err
N2fits  = [[[N2metal]],[[N2metal_err]],[[upper_tags]]] 

PA_spec    = sxpar(ha_hdr,'PA_SPEC')
scale      = sxpar(ha_hdr,'scale') ; in degrees
;scale      = scale*3600. ; in arcsecs


;=============================================
;2) Metallicity from Ha, NII, Hb, and OIII
if type eq "O3N2" then begin
OIIIcube = readfits(OIII_file,OIII_hdr)

;--------------------------------------------
;Align pixels from 2 sets. Align OIII to Halpha

writefits,'test.fits',oiiicube[*,*,1],oiii_hdr
oiii=readfits('test.fits',oiii1d_hdr)

hbeta     = oiiicube[*,*,6]
oiii_err  = oiiicube[*,*,9]
hbeta_err = oiiicube[*,*,11] 
velocity = oiiicube[*,*,0]
vel_disp = oiiicube[*,*,2]
bins     = oiiicube[*,*,3]
oiiisig  = oiiicube[*,*,4]
redshift = oiiicube[*,*,5]
hbsig    = oiiicube[*,*,7]
velsig   = oiiicube[*,*,8]
dispsig  = oiiicube[*,*,10]


hastrom,oiii,oiii1d_hdr,new_oiii,new_oiii_hdr,ha_hdr,missing=1./0.
hastrom,hbeta,oiii1d_hdr,new_hbeta,newhdr,ha_hdr,missing=1./0.
hastrom,oiii_err,oiii1d_hdr,new_oiii_err,newhdr,ha_hdr,missing=1./0.
hastrom,hbeta_err,oiii1d_hdr,new_hbeta_err,newhdr,ha_hdr,missing=1./0.
hastrom,velocity,oiii1d_hdr,new_velocity,newhdr,ha_hdr,missing=1./0.
hastrom,vel_disp,oiii1d_hdr,new_vel_disp,newhdr,ha_hdr,missing=1./0.
hastrom,bins,oiii1d_hdr,new_bins,newhdr,ha_hdr,missing=1./0.
hastrom,oiiisig,oiii1d_hdr,new_oiiisig,newhdr,ha_hdr,missing=1./0.
hastrom,redshift,oiii1d_hdr,new_redshift,newhdr,ha_hdr,missing=1./0.
hastrom,hbsig,oiii1d_hdr,new_hbsig,newhdr,ha_hdr,missing=1./0.
hastrom,velsig,oiii1d_hdr,new_velsig,newhdr,ha_hdr,missing=1./0.
hastrom,dispsig,oiii1d_hdr,new_dispsig,newhdr,ha_hdr,missing=1./0.

;Now all the new_param is aligned with the Halpha image
new_OIIIcube = [[[new_velocity]],[[new_oiii]],[[new_vel_disp]],[[new_bins]],[[new_oiiisig]],[[new_redshift]],[[new_hbeta]],[[new_hbsig]],[[new_velsig]],[[new_oiii_err]],[[new_dispsig]],[[new_hbeta_err]]]

namenewoiii = strmid(oiii_File,0,strpos(oiii_file,'.fits'))+'_aligned.fits'
writefits,namenewoiii,new_OIIIcube,new_oiii_hdr
;stop


;--------------------------------------------
;CALCULATE METALICITY
if where(new_OIII lt 0. and NII lt 0.) ne -1 then stop

O3N2 = alog10((new_OIII/new_hbeta)/(NII/halpha))
O3N2metal = 8.73-0.32*O3N2

O3N2_err = 0.43429*sqrt((new_OIII_err/new_OIII)^2+(new_hbeta_err/new_hbeta)^2+(NII_err/NII)^2+(halpha_err/halpha)^2)
O3N2metal_err = 0.32*O3N2_err
 
;; difference in metallicity using 2 methods
Diff_metal = O3N2metal-N2metal
Diff_metal_err = sqrt(O3N2metal_err^2+N2metal_err^2)

;;Halpha/Hbeta map
Ha_over_Hb = Halpha/new_Hbeta
Ha_over_Hb_err = sqrt((Halpha_err/Halpha)^2+(new_Hbeta_err/new_Hbeta)^2)


; MAKING BPT DIAGRAMS

N2divHa =  NII/halpha
O3divHb =  new_OIII/new_hbeta
N2divHa_err = sqrt((NII_err/NII)^2+(halpha_err/halpha)^2)*N2divHa
O3divHb_err = sqrt((new_OIII_err/new_OIII)^2+(new_hbeta_err/new_hbeta)^2)*O3divHb
pix_indices = indgen((size(nii))[1]*(size(nii))[2])

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
   map_agn = fltarr((size(nii))[1],(size(nii))[2])
   map_agn(*,*) = 1./0.
   goto, nobptdata
endif

  ; Find points that are above Kewley AGN line
AGN = 0.61/(alog10(n2divha)-0.47) + 1.19
agnpointsinn2divha = (where(alog10(O3divhb) ge agn)) ;indices in n2divha

if agnpointsinn2divha ne [-1] then begin   
   AGNpoints = pix_indices(where(alog10(O3divhb) ge agn)) ;indices in metal map
   map_agn = fltarr((size(nii))[1],(size(nii))[2])
   map_agn(AGNpoints) = 1.
   agnflag= 1
endif else begin 
   print, 'NO AGN. Yay!'
   map_agn = fltarr((size(nii))[1],(size(nii))[2])
   agnflag = 0
endelse

; Analyze indicated AGN points and make intensity vs wavelengths and fitted gaussian profiles
if agnflag eq 2 then begin
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

if agnflag eq 2 then oploterror,N2divHa(goodagninn2divha),O3divHb(goodagninn2divha),N2divHa_err(goodagninn2divha),o3divHb_err(goodagninn2divha),psym=4,thick=5

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

if agnflag eq 2 then begin 
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
 
;data_quality, path=path,fits_file=[ha_file,nii_file,oiii_file,hb_file],[1,2,1,2],[x_shift_ha,x_shift_ha,x_shift_hb,x_shift_hb],[y_shift_ha,y_shift_ha,y_shift_hb,y_shift_hb],size_x,size_y,nameplot=[name+' Halpha',name+' NII',name+' OIII',name+' Hbeta'],namegal=name



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Creating fits array
O3N2_fits = [[[N2metal]],[[O3N2metal]],[[Diff_metal]],[[Ha_over_hb]],[[map_agn]],[[N2metal_err]],[[O3N2metal_err]],[[Diff_metal_err]],[[Ha_over_Hb_err]]]

;MAKE HEADER
o3n2_header = ha_hdr
hdr_para_names = ['frame1','frame2','frame3','frame4','frame5','frame6','frame7','frame8','frame9']
hdr_para_vals = [1/0,1/0,1/0,1/0,1/0,1/0,1/0,1/0,1/0]
comment1 = ['N2metal','O3N2metal','O3N2metal-N2metal','Ha/Hb','Detected AGN points','N2metal Error','O3N2metal error','O3N2metal-N2metal error','Ha/Hb error']

for ii=0,n_elements(hdr_para_names)-1 do begin 
   sxaddpar,o3n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment1(ii)
endfor


writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_metallicity.fits',O3N2_fits,o3n2_header

;;;;;;;;;;;;;;;;;;
;Make Metallicity as a function of radius
if name eq 'nCSWA128' then xy_Hacont_peak = [22,35] else $
if name eq 'nCSWA165' then xy_Hacont_peak = [25,37] else goto, skip

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
   n2_header = ha_hdr
   hdr_para_names = ['frame1','frame2']
   hdr_para_vals = [1/0,1/0]
   comment = ['N2metal','N2metal Error']

   for ii=0,n_elements(hdr_para_names)-1 do sxaddpar,n2_header,hdr_para_names(ii),hdr_para_vals(ii),comment(ii)

   writefits, '/scr2/nichal/workspace/output/metallicity/'+name+'_N2_metallicity.fits',N2fits,n2_header
endif

end



