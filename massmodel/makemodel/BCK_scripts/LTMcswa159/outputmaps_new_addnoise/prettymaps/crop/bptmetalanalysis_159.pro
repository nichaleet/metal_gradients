pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((gradient*x+Zc-8.90)/0.57)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((gradient*x+Zc-8.62)/0.36)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcPP04_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((8.73-(gradient*x+Zc))/0.32)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = 10.^((8.66-(gradient*x+Zc))/0.28)
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro bptmetalanalysis_159,name,path,hadetect,nii,nii_err,ha,ha_err,oiii,oiii_err,hb,hb_err,outparams,angle=angle,xc=xc,yc=yc,pixelscale=pixelscale,slitwidth=slitwidth,inc=inc,pos=pos,minhadetect=minhadetect,restrict_rmax=restrict_rmax,noweight=noweight
;minhadetect is used in pixelate method
outdir = '/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps_new_addnoise/prettymaps/crop/metal_result/'
;set plot for eps maps
set_plot, 'ps'
page_height = 15.
page_width = 15.
plot_left = 3.
plot_bottom = 2.
xsize = 10.
ysize = 10.

;Define greek letters
alphaLetter = '!9' + String("141B) + '!X'  
betaletter  = '!9' + String("142B) + '!X'

;change position angle and inclination into radians
inc = inc*!pi/180.
pos = pos*!pi/180.
;stop
;load data
nii     = readfits(path+nii)
nii_err = readfits(path+nii_err)
ha      = readfits(path+ha)
ha_err  = readfits(path+ha_err)
hadetect = readfits(path+hadetect)
;stop
if file_search(path+oiii) ne '' then begin
   N2onlyflag = 0
   oiii    = readfits(path+oiii)
   oiii_err= readfits(path+oiii_err)
   hb      = readfits(path+hb)
   hb_err  = readfits(path+hb_err)
endif else begin
   N2onlyflag = 1
   oiii    = fltarr(size(ha,/dimension))+1./0.
   oiii_err= oiii
   hb      = oiii
   hb_err  = oiii
endelse
;stop
;fix the bad nii
bad = where(abs(nii) gt .1 or finite(nii_err) eq 0.)
nii(bad)     = 1./0.
nii_err(bad) = 1./0.
ha(bad)      = 1./0.
ha_err(bad)  = 1./0.

;calculate the n2, o3hb, and o3n2 indices
n2       = nii/ha
n2_err   = abs(n2)*sqrt((nii_err/nii)^2+(ha_err/ha)^2)
o3hb     = oiii/hb
o3hb_err = abs(o3hb)*sqrt((oiii_err/oiii)^2+(hb_err/hb)^2)

;o3n2 is usually negative when n2 is negative. So that means o3n2 is huge and should not be thought as zero
negnii = where(nii lt 0.)
nii_temp = nii
nii_err_temp = nii_err
nii_temp(negnii) = 0.5*(nii_temp(negnii)+nii_err_temp(negnii)) ;half way between 0 and maximum
toohigh = where(nii_temp gt 0.0015)
if toohigh(0) ne -1 then nii_temp(toohigh) = 0.0015
nii_err_temp(negnii) = nii_temp(negnii) 
n2_temp = nii/ha
n2_err_temp =  abs(n2_temp)*sqrt((nii_err_temp/nii_temp)^2+(ha_err/ha)^2)
o3n2     = o3hb/n2_temp
o3n2_Err = abs(o3n2)*sqrt((o3hb_err/o3hb)^2+(n2_err_temp/n2_temp)^2)
o3n2_err(where(o3n2 lt 0.)) = 1./0.
o3n2(where(o3n2 lt 0.)) = 1./0.

;stop
;write N2 map and O3N2 map
writefits,outdir+name+'_N2map.fits',[[[N2]],[[N2_err]]]
if N2onlyflag eq 0 then begin
   writefits,outdir+name+'_O3N2map.fits',[[[O3N2]],[[O3N2_err]]]
;plot BPT diagram
   bptdiagram,name,outdir,ha,n2,n2_err,o3hb,h2_err,o3hb_err,page_width,page_height,plot_left,plot_bottom,xsize,ysize
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Metalicity;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;1 The PP04,Steidel N2 index ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;Take care of no detection N2 business 
mederr = median(N2_err)
goodneg = where(N2 lt 0. and N2+N2_err gt 0.)
if goodneg(0) ne -1 then begin
   N2_err(goodneg) = 0.5*(N2(goodneg)+N2_err(goodneg))
   N2(goodneg) = N2_err(goodneg)
endif
badneg = where(N2 lt 0. and N2+N2_err lt 0.)
if badneg(0) ne -1 then begin
   N2(badneg) = 0.
   N2_err(badneg) = mederr
endif
;

;;1.1) Find Gradient from all individual pixels binned in annuli
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(minhadetect) eq 0. then minhadetect = 5.
badhadetection = where(hadetect lt minhadetect)
annuli_calc,N2,N2_err,xc,yc,pos,inc,pixelscale,badhadetection,radius,radius_pixel,N2_radius_pixel,N2_radius_err_pixel,radius_annuli,N2_radius_annuli,N2_radius_err_annuli
bad = where(N2_Radius_annuli lt 0.005)
if bad(0) ne -1 then remove,bad,radius_annuli,N2_radius_annuli,N2_radius_err_annuli
radius_annuli_of_N2 = radius_annuli
writefits,outdir+name+'_radius.fits',radius


;;1.2) Find Gradient along a pseudo slit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Get N2 as a function of distance along a slit and make a slitmaps
sizeim = size(N2,/dimension)
width = sizeim(1)
n_slits = fix(1.4*width/slitwidth) ;number of subslits
if name eq 'cswa165' then n_slits = fix(2.*width/slitwidth)

;Get N2 as a function of distance and get a slit map that shows a slit.
N2distance_slitmap,name,outdir,n_slits,slitwidth,pixelscale,xc,yc,angle,N2,N2_err,N2_distance,N2_distance_err,distance,N2mapslit,midpoint,slit_indices,slit_size
;The result fits map is in outdir+name+'_N2map.fits' with N2,N2_err, and N2mapslit

;Get the N2 index as a function of radius along major axis (N2_radius)
N2radius,name,outdir,N2mapslit,N2,N2_err,midpoint,slit_indices,slit_size,N2_radius,N2_radius_err,radius

;------------------------------PLOT---------------------------
;plot N2 as a function of radius
psname=outdir+name+'_N2radius.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

if min(N2_radius) lt 1.e-3 then yrange = [-3,alog10(max(N2_radius))] else yrange = alog10([min(N2_radius),max(N2_radius)])

;plot individual pixels
loadct,0
xrange=[0,max(radius)+2.]
if name eq 'cswa165' then xrange=[0,max(radius)+.2]
if name eq 'cswa19_maingal' then yrange=[0.005,.2]
;ploterror,radius_pixel,alog10(N2_radius_pixel),abs(0.434*N2_radius_err_pixel/N2_radius_pixel),nohat=nohat,psym=4,errcolor=200,yrange=yrange,xrange=xrange,xmargin=[10,7],ystyle=8,xtitle='radius(kpc)',ytitle='[NII]/H'+alphaletter,type=0,font=0,/nodata
;cgplot,radius_pixel,N2_radius_pixel,psym=14,color=200,/overplot
cgplot,radius_pixel,alog10(N2_radius_pixel),psym=200,yrange=yrange,xrange=xrange,xmargin=[10,7],ystyle=8,xtitle='radius(kpc)',ytitle='[NII]/H'+alphaletter,font=0 

;plot annuli
loadct,2
oploterror,radius_annuli,alog10(N2_radius_annuli),abs(0.434*N2_radius_err_annuli/N2_radius_annuli),psym=7,errcolor=5
cgplot,radius_annuli,alog10(N2_radius_annuli),psym=7,color=5,/overplot

;plot N2 along a pseudo slit
oploterror,radius,alog10(N2_radius),abs(0.434*N2_radius_err/N2_radius),psym=7,errcolor=200
cgplot,radius,alog10(N2_radius),psym=7,color=200,/overplot

;save the radius values
radius_of_n2 = radius

;add the second y axis in (12+logO/H)
metalrange=8.90+0.57*!y.crange
Axis,yaxis=1,yrange=metalrange,ylog=0,ytitle='12+log(O/H)!LPP04',font=0,ystyle=1

;----------- Fit Functions ---------------------------------------
;Fit the PP04 function along pseudo slit
A = [0.,8.90+0.57*alog10(N2_radius(0))]  ;guess parameter(gradient and Z center)
if A(1) lt 7. then A(1) = 7.2
weight = 1./N2_radius_err^2
x=radius_of_n2
y=N2_radius
;remove the outliers
outlier_lim=[]

outlier_lim = 0.009
if n_Elements(outlier_lim) ne 0 then begin   
   outlier1 = where(y lt outlier_lim)
   if outlier1(0) ne -1 then begin
      cgplot,x(outlier1),y(outlier1),psym=15,color=200,/overplot
      remove,outlier1,x,y,weight
   endif
endif

if keyword_set(noweight) then weight=weight*0.+1.


N2fit_pp04 = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',chisq=chisq,/noderivative)
A_PP04 = A
sigmaA_pp04 = sigmaA
;Fit the Steidel function along pseudo slit 
N2fit_steidel = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel',/noderivative)
A_steidel = A
sigmaA_steidel = sigmaA


;Fit the PP04 function for the annuli ;;;;;;;;;;;;;;;;;;;;;;;;;;;
weight = 1./N2_radius_err_annuli^2
if keyword_set(noweight) then weight=weight*0.+1.
x=radius_annuli_of_n2
y=N2_radius_annuli
if n_elements(restrict_rmax) ne 0 then begin
   rmax=max(radius_of_n2)
   bad = where(x gt rmax)
   cgplot,radius_annuli_of_n2(bad),N2_radius_annuli(bad),psym=15,color=5,/overplot
   remove,bad,x,y,weight
endif
;remove outliers
if n_Elements(outlier1) ne 0 then begin   
   if outlier1(0) ne -1 then begin
      cgplot,x(outlier1),y(outlier1),psym=15,color=5,/overplot
      remove,outlier1,x,y,weight
   endif
endif
N2fit_annuli_pp04 = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04',/noderivative,chisq=chisq,status=status)
;stop
A_annuli_PP04 = A
sigmaA_annuli_pp04 = sigmaA
;Fit the Steidel function for the annuli ;;;;;;;;;;;;;;;;;;;;;;;;;;;
N2fit_annuli_steidel = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel',/noderivative,status=status)
A_annuli_steidel = A
sigmaA_annuli_steidel = sigmaA

;make the fitted plot
xfit = !x.crange
funcPP04,xfit,A_PP04,yfit,pder
oplot,xfit,yfit,color=200
funcSteidel,xfit,A_steidel,yfit,pder
oplot,xfit,yfit,color=200
funcPP04,xfit,A_annuli_PP04,yfit,pder
oplot,xfit,yfit,color=5
funcSteidel,xfit,A_annuli_steidel,yfit,pder
oplot,xfit,yfit,color=5

device,/close

outparams=[A_pp04,A_steidel,A_annuli,A_annuli_steidel]

;plot N2 as a function of distance 
loadct,2
set_plot,'ps'
psname=outdir+name+'_N2distance.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
if min(N2_distance) lt 1.e-3 then yrange = [1.e-2,max(N2_distance)] else yrange = [min(N2_distance),max(N2_distance)] 

ploterror,distance,N2_distance,N2_distance_err,xtitle='distance(kpc)',ytitle='[NII]/H'+alphaletter,psym=7,type=1,font=0,errcolor=200,yrange=yrange,hatlength=!D.X_VSIZE/50.
  
device,/close


end
