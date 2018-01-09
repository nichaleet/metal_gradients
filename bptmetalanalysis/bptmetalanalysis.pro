pro funcPP04,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = (gradient*x+Zc-8.90)/0.57
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = (gradient*x+Zc-8.62)/0.36
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcPP04_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = (8.73-(gradient*x+Zc))/0.32
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro funcSteidel_O3N2,x,A,f,pder

gradient = A[0]
Zc = A[1]
f = (8.66-(gradient*x+Zc))/0.28
pder = fltarr(n_elements(x),n_elements(A)) ; no value returned

end

pro bptmetalanalysis,name,path,hadetect,nii,nii_err,ha,ha_err,oiii,oiii_err,hb,hb_err,angle=angle,xc=xc,yc=yc,pixelscale=pixelscale,slitwidth=slitwidth,inc=inc,pos=pos,minhadetect=minhadetect,restrict_rmax=restrict_rmax,noweight=noweigh,detectlim=detectlim
;minhadetect is used in pixelate method
outdir = '/scr2/nichal/workspace/output/bptmetalanalysis_result/'
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

;load data
nii     = readfits(path+nii)
nii_err = readfits(path+nii_err)
ha      = readfits(path+ha)
ha_err  = readfits(path+ha_err)
hadetect = readfits(path+hadetect)
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

;fix the bad nii
;fix too high
bad = where(abs(nii) gt .1 or finite(nii_err) eq 0.)
if bad(0) ne -1 then begin
   nii(bad)     = 1./0.
   nii_err(bad) = 1./0.
   ha(bad)      = 1./0.
   ha_err(bad)  = 1./0.
endif

;fix good negative 
goodneg = where(nii lt 0. and nii+nii_err gt 0.)
if goodneg(0) ne -1 then begin
   nii(goodneg) = 0.5*(nii(goodneg)+nii_err(goodneg))
   nii_err(goodneg) = nii(goodneg)
endif

;fix bad negative
badneg = where(nii lt 0. and nii+nii_err lt 0.)
if badneg(0) ne -1 then begin
   nii(badneg) =  0.5*min(Nii(where(Nii gt 0.))) ;half of the smallest better detection
   nii_err(badneg) = detectlim
endif

;calculate the n2, o3hb, and o3n2 indices
n2       = nii/ha
n2_err   = abs(n2)*sqrt((nii_err/nii)^2+(ha_err/ha)^2)
o3hb     = oiii/hb
o3hb_err = abs(o3hb)*sqrt((oiii_err/oiii)^2+(hb_err/hb)^2)

o3n2     = o3hb/n2
o3n2_Err = abs(o3n2)*sqrt((o3hb_err/o3hb)^2+(n2_err/n2)^2)
o3n2_err(where(o3n2 le 0.)) = 1./0.
o3n2(where(o3n2 le 0.)) = 1./0.

logN2 = alog10(N2)
logN2(where(finite(logN2) eq 0)) = 1./0. ; fix the -inf to +inf
logN2_err = abs(0.4343*N2_err/N2)

logo3hb = alog10(o3hb)
logo3hb(where(finite(logo3hb) eq 0)) = 1./0.
logo3hb_err = abs(0.4343*o3hb_err/o3hb)

logO3N2 = alog10(O3N2)
logO3N2(where(finite(logO3N2) eq 0)) = 1./0.
logO3N2_err = (0.4343*O3N2_err/O3N2)

;write N2 map and O3N2 map
writefits,outdir+name+'_N2map.fits',[[[N2]],[[N2_err]]]
if N2onlyflag eq 0 then begin
   writefits,outdir+name+'_O3N2map.fits',[[[O3N2]],[[O3N2_err]]]
;plot BPT diagram
   bptdiagram,name,outdir,logn2,logn2_err,logo3hb,logo3hb_err,page_width,page_height,plot_left,plot_bottom,xsize,ysize
endif
stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Metalicity;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;1 The PP04,Steidel N2 index ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;1.1) Find Gradient from all individual pixels binned in annuli
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(minhadetect) eq 0. then minhadetect = 5.
badhadetection = where(hadetect lt minhadetect)
;stop
annuli_calc,logN2,logN2_err,xc,yc,pos,inc,pixelscale,badhadetection,radius,radius_pixel,N2_radius_pixel,N2_radius_err_pixel,radius_annuli,N2_radius_annuli,N2_radius_err_annuli
bad = where(N2_Radius_annuli lt -3.)
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
N2distance_slitmap,name,outdir,n_slits,slitwidth,pixelscale,xc,yc,angle,logN2,logN2_err,N2_distance,N2_distance_err,distance,N2mapslit,midpoint,slit_indices,slit_size
;The result fits map is in outdir+name+'_N2map.fits' with N2,N2_err, and N2mapslit

;Get the N2 index as a function of radius along major axis (N2_radius)
N2radius,name,outdir,N2mapslit,logN2,logN2_err,midpoint,slit_indices,slit_size,N2_radius,N2_radius_err,radius

;------------------------------PLOT---------------------------
;plot N2 as a function of radius
psname=outdir+name+'_N2radius.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

if min(N2_radius) lt -6 then yrange = [-3,max(N2_radius)] else yrange = [min(N2_radius),max(N2_radius)]
yrange = [-2.3,-0.5]
if name eq 'cswa31' then yrange = [-1.5,0]
if name eq 'cswa128' or name eq 'cswa159' then yrange = [-2.,0]
if name eq 'cswa165' then yrange=[-1.5,0.]
;plot individual pixels
loadct,0
if keyword_set(restrict_rmax) then begin
   xrange = [0,restrict_rmax] 
   if (where(radius_pixel gt restrict_rmax))[0] ne -1 then remove,where(radius_pixel gt restrict_rmax),radius_pixel,n2_radius_pixel,n2_radius_err_pixel
   if (where(radius_annuli_of_n2 gt restrict_rmax))[0] ne -1 then remove, where(radius_annuli_of_n2 gt restrict_rmax),radius_annuli_of_n2,n2_radius_annuli,n2_radius_err_annuli
   if (where(radius gt restrict_rmax))[0] ne -1 then remove, where(radius gt restrict_rmax),radius,n2_radius,n2_radius_err
endif else xrange=[0,max(radius)+2.]
if name eq 'cswa28' then xrange = [0,8.]
if name eq 'cswa165' then xrange=[0,1.5]
if name eq '1148' then xrange=[0,2.]
if name eq '0744' then xrange=[0,1.2]
if name eq '0744' or name eq '1148' then yrange=[-2.,0.]

ploterror,radius_pixel,N2_radius_pixel,N2_radius_err_pixel,nohat=nohat,psym=4,errcolor=255,yrange=yrange,xrange=xrange,xmargin=[10,7],ystyle=8,xtitle='radius(kpc)',ytitle='[NII]/H'+alphaletter,type=0,font=0,/nodata
cgplot,radius_pixel,N2_radius_pixel,psym=14,color=200,/overplot

;plot annuli
loadct,2
oploterror,radius_annuli_of_n2,N2_radius_annuli,N2_radius_err_annuli,psym=3,errcolor=5
cgplot,radius_annuli,N2_radius_annuli,psym=20,color=5,/overplot

;plot N2 along a pseudo slit
if name eq 'cswa19_maingal' and (where(radius gt 1.8))[0] ne -1 then remove,where(radius gt 1.8), radius,N2_radius,N2_radius_err

oploterror,radius,N2_radius,N2_radius_err,psym=3,errcolor=200
cgplot,radius,N2_radius,psym=16,color=200,/overplot

;save the radius values
radius_of_n2 = radius

;add the second y axis in (12+logO/H)
metalrange=8.90+0.57*!y.crange
Axis,yaxis=1,yrange=metalrange,ytitle='12+log(O/H)!LPP04',font=0,ystyle=1;,ylog=0


;----------- Fit Functions ---------------------------------------
;Fit the PP04 function along pseudo slit
A = [0.,8.90+0.57*N2_radius(0)]  ;guess parameter(gradient and Z center)
if A(1) lt 7. then A(1) = 7.2
weight = 1./N2_radius_err^2
x=radius_of_n2
y=N2_radius

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
;oplot,xfit,yfit,color=200
funcPP04,xfit,A_annuli_PP04,yfit,pder
oplot,xfit,yfit,color=5
funcSteidel,xfit,A_annuli_steidel,yfit,pder
;oplot,xfit,yfit,color=5

device,/close

;Append the metal gradient values to the N2method gradient file-------------------
openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, ''
printf,1, name
printf,1, systime()
printf,1, 'N2 PP04 Pseudo slit'
printf,1, 'metal gradient =',A_pp04(0),'+/-',sigmaA_pp04(0)
printf,1, 'central metal =', A_pp04(1),'+/-',sigmaA_pp04(1)
printf,1, 'N2 Steidel Pseudo slit'
printf,1, 'metal gradient =',A_steidel(0),'+/-',sigmaA_steidel(0)
printf,1, 'central metal =', A_steidel(1),'+/-',sigmaA_steidel(1)
printf,1, 'N2 PP04 Annuli'
printf,1, 'metal gradient =',A_annuli_pp04(0),'+/-',sigmaA_annuli_pp04(0)
printf,1, 'central metal =', A_annuli_pp04(1),'+/-',sigmaA_annuli_pp04(1)
printf,1, 'N2 Steidel Annuli'
printf,1, 'metal gradient =',A_annuli_steidel(0),'+/-',sigmaA_annuli_steidel(0)
printf,1, 'central metal =', A_annuli_steidel(1),'+/-',sigmaA_annuli_steidel(1)
close,1

;plot N2 as a function of distance 
loadct,2
set_plot,'ps'
psname=outdir+name+'_N2distance.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
if min(N2_distance) lt -3 then yrange = [-2,max(N2_distance)] else yrange = [min(N2_distance),max(N2_distance)] 

ploterror,distance,N2_distance,N2_distance_err,xtitle='distance(kpc)',ytitle='[NII]/H'+alphaletter,psym=7,type=0,font=0,errcolor=200,yrange=yrange,hatlength=!D.X_VSIZE/50.
  
device,/close

if N2onlyflag eq 0 then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;2 The PP04&Steidel O3N2 index ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;2.1) Find Gradient from all individual pixels of O3N2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

annuli_calc,logO3N2,logO3N2_err,xc,yc,pos,inc,pixelscale,badhadetection,radius,radius_pixel,O3N2_radius_pixel,O3N2_radius_err_pixel,radius_annuli,O3N2_radius_annuli,O3N2_radius_err_annuli

;;2.2) Find Gradient along a pseudo slit O3N2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Get O3N2 as a function of distance and radius along a slit and make a slitmaps

o3n2_distance_radius,logO3N2,logO3N2_err,n2mapslit,n_slits,midpoint,slit_size,o3n2mapslit,distance,o3n2_distance,o3n2_distance_err,radius,o3n2_radius,o3n2_radius_err

;write O3N2 map 
writefits,outdir+name+'_O3N2map.fits',[[[logO3N2]],[[logO3N2_err]],[[O3N2mapslit]]]

;------------------------------PLOT---------------------------

;plot O3N2 as a function of distance -----------------------
loadct,2
set_plot,'ps'
psname=outdir+name+'_O3N2distance.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   
ploterror,distance,O3N2_distance,O3N2_distance_err[*],xtitle='distance(kpc)',ytitle='([OIII]/H'+betaletter+')/([NII]/H'+alphaletter+')',psym=7,type=0,font=0,errcolor=200,yrange=[1.,3.],HATLENGTH=!D.X_VSIZE/50.
   
device,/close

;plot O3N2 as a function of radius-----------------------
psname=outdir+name+'_O3N2radius.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

if max(O3N2_radius) gt 2.5 then yrange = [min(O3N2_radius),3.] else yrange = [min(O3N2_radius),max(O3N2_radius)] 
yrange = [3.5,1.]
if name eq 'cswa128' or name eq 'cswa159' then yrange = [3.,0.5]
if name eq 'cswa165' then yrange = [2.5,0.]
;plot individual pixels
loadct,0
if keyword_set(restrict_rmax) then begin
   xrange=[0,restrict_rmax] 
   if (where(radius_pixel gt restrict_rmax))[0] ne -1 then remove,where(radius_pixel gt restrict_rmax),radius_pixel,o3n2_radius_pixel,o3n2_radius_err_pixel
   if (where(radius_annuli gt restrict_rmax))[0] ne -1 then remove, where(radius_annuli gt restrict_rmax),radius_annuli,o3n2_radius_annuli,o3n2_radius_err_annuli
   if (where(radius gt restrict_rmax))[0] ne -1 then remove, where(radius gt restrict_rmax),radius,o3n2_radius,o3n2_radius_err
endif else xrange=[0,max(radius)+2.]

if name eq 'cswa165' then xrange=[0,1.5]
ploterror,radius_pixel,O3N2_radius_pixel,O3N2_radius_err_pixel,nohat=nohat,psym=4,errcolor=255,yrange=yrange,xrange=xrange,xmargin=[10,7],ystyle=8,xtitle='radius(kpc)',ytitle='([OIII]/H'+betaletter+')/([NII]/H'+alphaletter+')',type=0,font=0,/nodata
cgplot,radius_pixel,O3N2_radius_pixel,psym=14,color=200,/overplot

;plot annuli
loadct,2
oploterror,radius_annuli,O3N2_radius_annuli,O3N2_radius_err_annuli[*],psym=3,errcolor=5
cgplot,radius_annuli,O3N2_radius_annuli,psym=20,color=5,/overplot

;plot O3N2 along a pseudo slit
oploterror,radius,O3N2_radius,O3N2_radius_err,psym=3,errcolor=200
cgplot,radius,O3N2_radius,psym=16,color=200,/overplot

;add the second y axis in (12+logO/H)
metalrange=8.73-0.32*!y.crange
Axis,yaxis=1,yrange=metalrange,ylog=0,ytitle='12+log(O/H)!L PP04',font=0,ystyle=1
   
;save the radius values
radius_of_o3n2 = radius

;----------- Fit Functions ---------------------------------------
;Fit the PP04 function along pseudo slit
   A = [0.,8.66-0.28*O3N2_radius(0)];guess parameter(gradient and Z center)
if n_elements(o3n2_radius) gt 2. then begin
   weight = 1./O3N2_radius_err^2
   ;if name eq 'cswa165' then weight = fltarr(n_elements(weight))+1.
   x = radius_of_o3n2
   y = O3N2_radius
   O3N2fit_pp04 = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04_O3N2',/noderivative)
   A_PP04 = A
   sigmaA_pp04 = sigmaA

;Fit the Steidel function along pseudo slit
   A = [0.,8.66-0.28*O3N2_radius(0)];guess parameter(gradient and Z center)
   O3N2fit_steidel = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel_O3N2',/noderivative)
   A_steidel = A
   sigmaA_steidel = sigmaA

;make the fitted plot
   xfit = !x.crange
   funcPP04_O3N2,xfit,A_PP04,yfit,pder
   oplot,xfit,yfit,color=200
   funcSteidel_O3N2,xfit,A_steidel,yfit,pder
;   oplot,xfit,yfit,color=200
endif else begin
   A_PP04=[0,0]
   sigmaA_pp04 = [0,0]
   A_steidel = [0,0]
   sigmaA_steidel =[0,0]
endelse

;Fit the PP04 function for the annuli 
   weight = 1./O3N2_radius_err_annuli^2
   if keyword_set(noweight) then weight=weight*0.+1.
   ;if name eq 'cswa165' then weight = fltarr(n_elements(weight))+1.
   x=radius_annuli
   y=O3N2_radius_annuli

   O3N2fit_annuli_pp04 = curvefit(x,y,weight,A,sigmaA,function_name='funcPP04_O3N2',/noderivative)
   A_annuli_PP04 = A
   sigmaA_annuli_pp04 = sigmaA
;Fit the Steidel function for the annuli

   O3N2fit_annuli_steidel = curvefit(x,y,weight,A,sigmaA,function_name='funcSteidel_O3N2',/noderivative)
   A_annuli_steidel = A
   sigmaA_annuli_steidel = sigmaA

;make the fitted plot
   xfit = !x.crange
   funcPP04_O3N2,xfit,A_annuli_PP04,yfit,pder
   oplot,xfit,yfit,color=5
   funcSteidel_O3N2,xfit,A_annuli_steidel,yfit,pder
;   oplot,xfit,yfit,color=5
   device,/close
stop
;Append the metal gradient values to the N2method gradient file-------------------
   openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
   printf,1, ''
   printf,1, 'O3N2 PP04 Pseudo slit'
   printf,1, 'metal gradient =',A_pp04(0),'+/-',sigmaA_pp04(0)
   printf,1, 'central metal =', A_pp04(1),'+/-',sigmaA_pp04(1)
   printf,1, 'O3N2 Steidel Pseudo slit'
   printf,1, 'metal gradient =',A_steidel(0),'+/-',sigmaA_steidel(0)
   printf,1, 'central metal =', A_steidel(1),'+/-',sigmaA_steidel(1)
   printf,1, 'O3N2 PP04 Annuli'
   printf,1, 'metal gradient =',A_annuli_pp04(0),'+/-',sigmaA_annuli_pp04(0)
   printf,1, 'central metal =', A_annuli_pp04(1),'+/-',sigmaA_annuli_pp04(1)
   printf,1, 'O3N2 Steidel Annuli'
   printf,1, 'metal gradient =',A_annuli_steidel(0),'+/-',sigmaA_annuli_steidel(0)
   printf,1, 'central metal =', A_annuli_steidel(1),'+/-',sigmaA_annuli_steidel(1)
   close,1
endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;3 N2 Maiolino08: ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;3.1) Individual pixel
;------------------------------------------
;Calculate metallicity at each pixel if the fits file doesn't exist.
m08exist = file_search('/scr2/nichal/workspace/output/bptmetalanalysis_result/',+name+'_m08_metal_log.fits')
if m08exist eq '' then begin
   set_plot,'x'
   Bayesian_metal_log,name,logo3hb,logo3hb_err,logn2,logn2_err,m08metal,median_m08metal,lowerbound_m08metal,upperbound_m08metal,typebayesian
   m08metal_err = (upperbound_m08metal-lowerbound_m08metal)/2.
   writefits,'/scr2/nichal/workspace/output/bptmetalanalysis_result/'+name+'_m08_metal_log.fits',[[[m08metal]],[[median_m08metal]],[[lowerbound_m08metal]],[[upperbound_m08metal]],[[m08metal_err]],[[typebayesian]]]
endif else begin
   m08 = readfits('/scr2/nichal/workspace/output/bptmetalanalysis_result/'+name+'_m08_metal_log.fits')
   m08metal            = m08[*,*,0]
   median_m08metal     = m08[*,*,1]
   lowerbound_m08metal = m08[*,*,2]
   upperbound_m08metal = m08[*,*,3]
   m08metal_err        = m08[*,*,4]
   typebayesian        = m08[*,*,5]
endelse
;------------------------------------------
annuli_calc,m08metal,m08metal_err,xc,yc,pos,inc,pixelscale,badhadetection,radius,radius_pixel,m08metal_radius_pixel,m08metal_radius_err_pixel,radius_annuli,m08metal_radius_annuli,m08metal_radius_err_annuli

;3.2) N2 along pseudo slit
Bayesian_metal_slit_n2,name,n2_radius,n2_radius_err,radius_of_n2,m08_n2_metal,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_n2_metal_err

;3.3) N2 along Annuli
Bayesian_metal_slit_n2,name,n2_radius_annuli,n2_radius_err_annuli,radius_annuli_of_n2,m08_n2_metal_annuli,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_n2_metal_err_annuli

;------------------------------PLOT---------------------------

;plot m08 as a function of radius
psname=outdir+name+'_m08_n2.eps'
device, filename = psname,xsize = 10,ysize = 8, $
xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
yrange= [min([m08_N2_metal_annuli,m08metal_radius_pixel]),max([m08_N2_metal_annuli,m08metal_radius_pixel])]
;plot individual pixels
loadct,0
xrange=[0,max(radius_pixel)+2.]
if name eq 'cswa165' then xrange=[0,max(radius_pixel)+.2]
ploterror,radius_pixel,m08metal_radius_pixel,m08metal_radius_err_pixel,nohat=nohat,psym=4,errcolor=200,yrange=yrange,xrange=xrange,xmargin=[10,7],xtitle='radius(kpc)',ytitle='12+log(O/H)!LM08,N2',font=0,/nodata
cgplot,radius_pixel,m08metal_radius_pixel,psym=14,color=200,/overplot

;plot annuli
loadct,2
oploterror,radius_annuli_of_n2,m08_n2_metal_annuli,m08_n2_metal_err_annuli,psym=7,errcolor=5
cgplot,radius_annuli_of_n2,m08_n2_metal_annuli,psym=7,color=5,/overplot

;plot pseudo slit
oploterror,radius_of_n2,m08_n2_metal,m08_n2_metal_err,psym=7,errcolor=200
cgplot,radius_of_n2,m08_n2_metal,psym=7,color=200,/overplot

;----------- Fit Functions ---------------------------------------
;Fit the function along pseudo slit
A_slit=linfit(radius_of_n2,m08_n2_metal,measure_errors=m08_n2_metal_err,sigma=sigmaA_slit)
;Fit the function along annuli
A_annuli=linfit(radius_annuli_of_n2,m08_n2_metal_annuli,measure_errors=m08_n2_metal_err_annuli,sigma=sigmaA_annuli)

;make the fitted plot
xfit = !x.crange
yfit = A_annuli(0)+A_annuli(1)*xfit
oplot,xfit,yfit,color=5
yfit = A_slit(0)+A_slit(1)*xfit
oplot,xfit,yfit,color=200
device,/close
;Append the metal gradient values to the gradient file

openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
printf,1, ''
printf,1, 'N2 M08 Slit'
printf,1, 'metal gradient =',A_slit(1),'+/-',sigmaA_slit(1)
printf,1, 'central metal =', A_slit(0),'+/-',sigmaA_slit(0)
printf,1, 'N2 M08 Annuli'
printf,1, 'metal gradient =',A_annuli(1),'+/-',sigmaA_annuli(1)
printf,1, 'central metal =', A_annuli(0),'+/-',sigmaA_annuli(0)
close,1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;4 O3N2 Maiolino08: ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if N2onlyflag eq 0 then begin
;4.1) o3hb and o3n2_M08 metallicity along pseudo slits
   o3hb_distance_radius,logo3hb,logo3hb_err,n2mapslit,n_slits,midpoint,slit_size,o3hbmapslit,distance,o3hb_distance,o3hb_distance_err,radius,o3hb_radius,o3hb_radius_err
   radius_of_o3hb=radius
   distance_of_o3hb = distance
;write o3hb map 
   writefits,outdir+name+'_o3hbmap.fits',[[[o3hb]],[[o3hb_err]],[[o3hbmapslit]]]
   
;find the match of radius
   match,radius_of_n2,radius_of_o3hb,sub_n2,sub_o3hb
   radius_match = radius_of_n2(sub_n2)
   o3hb_radius_match = o3hb_radius(sub_o3hb)
   o3hb_radius_err_match = o3hb_radius_err(sub_o3hb)
   n2_radius_match   = n2_radius(sub_n2)
   n2_radius_err_match   = n2_radius_err(sub_n2)
   
;calculate the metalicity
   Bayesian_metal_slit_o3n2,name,n2_radius_match,n2_radius_err_match,o3hb_radius_match,o3hb_radius_err_match,radius_match,m08_o3n2_metal,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_o3n2_metal_err
 
;4.2) o3hb and o3n2_M08 metallicity annuli
   annuli_calc,logO3hb,logO3hb_err,xc,yc,pos,inc,pixelscale,badhadetection,radius,radius_pixel,o3hb_radius_pixel,o3hb_radius_err_pixel,radius_annuli_of_o3hb,o3hb_radius_annuli,o3hb_radius_err_annuli
   testbad=where(finite(o3hb_radius_annuli) eq 0.)
   if testbad(0) ne -1 then remove,testbad,radius_annuli_of_o3hb,o3hb_radius_annuli,o3hb_radius_err_annuli
;find the match of radius
   match,radius_annuli_of_n2,radius_annuli_of_o3hb,sub_n2,sub_o3hb
   radius_annuli_match = radius_annuli_of_n2(sub_n2)
   o3hb_radius_annuli_match = o3hb_radius_annuli(sub_o3hb)
   o3hb_radius_err_annuli_match = o3hb_radius_err_annuli(sub_o3hb)
   n2_radius_annuli_match   = n2_radius_annuli(sub_n2)
   n2_radius_err_annuli_match   = n2_radius_err_annuli(sub_n2)
;calculate the metalicity
   Bayesian_metal_slit_o3n2,name,n2_radius_annuli_match,n2_radius_err_annuli_match,o3hb_radius_annuli_match,o3hb_radius_err_annuli_match,radius_annuli_match,m08_o3n2_metal_annuli,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,m08_o3n2_metal_err_annuli

;------------------------------PLOT---------------------------
;-----PLOT O3HB-----
;plot o3hb as a function of distance 
   loadct,2
   set_plot,'ps'
   psname=outdir+name+'_o3hbdistance.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   
   ploterror,distance_of_o3hb,o3hb_distance,o3hb_distance_err[*],xtitle='distance(kpc)',ytitle='[OIII]/H'+betaletter,psym=7,type=0,font=0,errcolor=200,yrange=[min(o3hb_distance),max(o3hb_distance)],HATLENGTH=!D.X_VSIZE/50.
   
   device,/close

;plot o3hb as a function of radius
   psname=outdir+name+'_o3hbradius.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   ploterror,radius_of_o3hb,o3hb_radius,o3hb_radius_err[*],xtitle='radius(kpc)',ytitle='log([OIII]/H'+betaletter+')',psym=7,type=0,font=0,errcolor=200,yrange=[min(o3hb_radius),max(o3hb_radius)],xmargin=[10,7],HATLENGTH=!D.X_VSIZE/50.
   device,/close   

;-----PLOT O3N2-----
;plot m08 as a function of radius
   psname=outdir+name+'_m08_o3n2.eps'
   device, filename = psname,xsize = 10,ysize = 8, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

;plot individual pixels o3n2 metallicity Maiolino
   loadct,0
   ploterror,radius_pixel,m08metal_radius_pixel,m08metal_radius_err_pixel,nohat=nohat,psym=4,errcolor=200,xrange=[0,max(radius_of_n2)+.2],xmargin=[10,7],xtitle='radius(kpc)',ytitle='12+log(O/H)!LM08,O3N2',font=0,/nodata
   cgplot,radius_pixel,m08metal_radius_pixel,psym=14,color=200,/overplot

;plot annuli o3n2 metallicity Maiolino
   loadct,2
   oploterror,radius_annuli_match,m08_o3n2_metal_annuli,m08_o3n2_metal_err_annuli,psym=7,errcolor=5
   cgplot,radius_annuli_match,m08_o3n2_metal_annuli,psym=7,color=5,/overplot
   
;plot pseudo slit
   oploterror,radius_match,m08_o3n2_metal,m08_o3n2_metal_err,psym=7,errcolor=200
   cgplot,radius_match,m08_o3n2_metal,psym=7,color=200,/overplot

;-------------------------------------------------------

;---------------Fit the function---------------

;fit the O3N2 M08 along pseudo slit
   A_slit=linfit(radius_match,m08_o3n2_metal,measure_errors=m08_o3n2_metal_err,sigma=sigmaA_slit)
;fit the O3N2 M08 along annuli
   A_annuli = linfit(radius_annuli_match,m08_o3n2_metal_annuli,measure_errors=m08_o3n2_metal_err_annuli,sigma=sigmaA_annuli)

;make the fitted plot
   xfit = !x.crange
   yfit = A_annuli(0)+A_annuli(1)*xfit
   oplot,xfit,yfit,color=5
   yfit = A_slit(0)+A_slit(1)*xfit
   oplot,xfit,yfit,color=200
   device,/close

;Append the metal gradient values to the gradient file
   openw, 1,'/scr2/nichal/workspace/output/bptmetalanalysis_result/gradient.dat',/append
   printf,1, 'O3N2 M08 Slit'
   printf,1, 'metal gradient =',A_slit(1),'+/-',sigmaA_slit(1)
   printf,1, 'central metal =', A_slit(0),'+/-',sigmaA_slit(0)
   printf,1, 'O3N2 M08 Annuli'
   printf,1, 'metal gradient =',A_annuli(1),'+/-',sigmaA_annuli(1)
   printf,1, 'central metal =', A_annuli(0),'+/-',sigmaA_annuli(0)
   close,1
endif
stop
end
