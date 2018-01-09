pro makerotmap,name,angle,xref,yref,xmid,ymid,x,y,n_slits,velocityfield,velocityfielderr,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,veldisptag,slit_indices
  slitmap = velocityfield ;this map will show where the slits are major axis
  ;equation for major axis is y=mx+c st. m=slope_major
  m_major = tan(angle/180.*!pi)  
  c_major = ymid-m_major*xmid
  ;distance of a point(x,y) to this major axis is |-mx+y-c|/sqrt(m^2+1)
  dist = abs((-1.*m_major*x)+y-c_major)/sqrt(m_major^2+1.)
  
  ;Next are the angle and equations for the slit dividers 
  angle_perpend = angle-90.  ;angle of the little slit
  slope = tan(angle_perpend/180.*!pi)  

; y = slope*x+intercept
; intercept = y-slope*x
  intercept = yref-slope(0)*xref

  velocity = fltarr(n_slits)
  velocity_err = fltarr(n_slits)
  slit_indices = findgen(n_slits)
;  xplot = fltarr(n_slits)
  for i=0,n_slits-1 do begin
     distmax = 1.5
     if name eq '0744' or name eq '1148' then distmax = 3.
     pix_in_slit = where(y le slope*x+intercept(i+1) and y ge slope*x+intercept(i) and dist le distmax) ; The slit width is ~3 pixels  
     if pix_in_slit(0) ne -1 then begin
        slitmap(pix_in_slit)=slit_indices(i)*100.
;        yslit = median(y(pix_in_slit),/even)
;        xslit = median(x(pix_in_slit),/even)
;        xplot = sqrt((xslit-xmid)^2+(yslit-ymid)^2)
        goodvel = velocityfield(pix_in_slit) 
        goodvelerr = velocityfielderr(pix_in_slit)
        thegood = where(finite(goodvel) and finite(goodvelerr))
        print,'There are', n_Elements(thegood),'pixels with velocity in this slit'
        simplemean = mean(goodvel(thegood))
        meanerr,goodvel(thegood),goodvelerr(thegood),wmean,sigmam,sigmad
        print,simplemean,wmean,' simple mean and weighted mean'

        if n_elements(thegood) eq 1 then begin
           velocity(i) = simplemean
           velocity_err(i) = stddev(goodvel(thegood))
        endif else begin
           velocity(i) = wmean
           velocity_err(i) = sigmad 
        endelse

     endif else print,'There is no pixel in this slit.'
  endfor
  stop
  if veldisptag ne 1 then writefits,'/scr2/nichal/workspace/slitmaps/slitmap.fits',[[[slitmap]],[[velocityfield]]]
  xplot = findgen(n_slits)*slit_size

  badpix=where(finite(velocity_err) eq 0 or velocity_err eq 0.) ;This will remove all the subslit with 1 pixel in it.
  if name eq 'cswa11_secondgal' then badpix = where(finite(velocity_err) eq 0 or velocity_err eq 0. or velocity_err gt 80. or velocity lt -170.)
  if badpix(0) ne -1 then remove,badpix, xplot,velocity,velocity_err,slit_indices

  yplot=velocity
  yploterr=velocity_err
  
  ploterror,xplot,velocity,velocity_err,xtitle='distance(kpc)',ytitle='velocity (km/s)',psym=1,title='angle='+string(angle)
  currentmaxvel = max(velocity(where(finite(velocity) eq 1)))

  rotationvelocity = currentmaxvel-min(velocity)
  velerr1 = velocity_err(where(velocity eq currentmaxvel))
  velerr2 = velocity_err(where(velocity eq min(velocity)))
  rotationvelocityerr = sqrt(velerr1^2+velerr2^2)
  rotationvelocityerr = rotationvelocityerr(0)
  ;wait,.5
end

;=======================================================================



;=======================================================================
;=======================================================================
;=======================================================================
pro rotationcurve_fix,path,velfile,velerrfile,dispfile,disperrfile,hadetectionfile,xc,yc,bestguess_angle,pixelscale,N2mapfile,N2errmapfile,M08mapfile,M08errmapfile,name
;Mar2015
;Input files: velocity map and NII/Ha map
;example for cswa128
;rotationcurve,'/scr2/nichal/workspace/output/cropvel/cswa128.fits',14,14,135.,0.04,2.22,'/scr2/nichal/workspace/output/cropmetal/cswa128.fits'
;pixelscale is in kpc per pixel
;bestguess_angle is best guest angle for major axis. It's the angle to x-axis counter clockwise

hadetection = readfits(path+hadetectionfile)
velocityfield= readfits(path+velfile)
velocityfielderr = readfits(path+velerrfile)
dispmap = readfits(path+dispfile)
dispmaperr=readfits(path+disperrfile)
N2map = readfits(path+N2mapfile,hdr)
N2maperr=readfits(path+N2errmapfile)

below5 = where(hadetection le 5.)
velocityfield(below5) = 1./0.
velocityfielderr(Below5)= 1./0.
dispmap(below5) = 1./0.
dispmaperr(below5) = 1./0.
N2map(below5) = 1./0.
N2maperr(below5) = 1./0.

sizeim = size(velocityfield)
width = sizeim(2)
read,slitwidth,prompt='slitwidth: '
;slit width in pixel -> This is the size of the subslit that divide the radius into radius bin, not the size of the big slit.
;*********** Change slit width here!!!

n_slits = fix(1.4*width/slitwidth)
if name eq 'cswa165' then n_slits = fix(2.*width/slitwidth)
;Calculate physical size per pixel
;theta_1kpc = zang(1.,z) ; in arcsecond
;pix_size  = pixelscale/theta_1kpc ;in kpc per pixel
pix_size = pixelscale ; in kpc per pixel


xmid = xc
ymid = yc

;xref, yref are the points that each slits are defined 
if bestguess_angle ge 45. and bestguess_angle le 135. then begin
   xref = fltarr(n_slits+1)+xmid    ; xref = xmid
   yref = findgen(n_slits+1)*slitwidth ;yref = 0,slitwidt,2slitwidth,...
   slit_size = slitwidth*pix_size*sin(bestguess_angle*!pi/180.)
endif
if bestguess_angle gt 135. then begin
   yref = fltarr(n_slits+1)+ymid       ; yref = ymid
   xref = reverse(findgen(n_slits+1)*slitwidth) ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
   slit_size = -1.*slitwidth*pix_size*cos(bestguess_angle*!pi/180.)
endif
if bestguess_angle lt 45. then begin
   yref = fltarr(n_slits+1)+ymid                ; yref = ymid
   xref = findgen(n_slits+1)*slitwidth ; xref = n*slitwidth,(n-1)*slitwidth,...,slitwidth,0
   slit_size = slitwidth*pix_size*cos(bestguess_angle*!pi/180.)
endif

if slit_size lt 0. then stop,'Stop: Slit size is negative'

;stop
ind = array_indices(velocityfield,findgen(n_elements(velocityfield)-1))
x = ind(0,*)
y = ind(1,*)


set_plot,'x'

makerotmap,name,bestguess_angle,xref,yref,xmid,ymid,x,y,n_slits,velocityfield,velocityfielderr,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,0,slit_indices
maxvel = rotationvelocity
maxvelerr = rotationvelocityerr

print, 'Max rotation velosity  is',maxvel,'+/-',maxvelerr,'with angle',bestguess_angle,'degree'


file_mkdir,'/scr2/nichal/workspace/output/rotationcurve1D_fromKineticmodel'
pos1 = strpos(velfile,'/',/reverse_Search)
pos2 = strpos(velfile,'.fits')
loadct,39

;shift the median velocity to 0

slitmap = readfits('/scr2/nichal/workspace/slitmaps/slitmap.fits')
midindex= slitmap(xc,yc,0)/100.
print, 'Midindex is ',midindex
midpoint=where(slit_indices eq midindex)
midpoint = midpoint(0)
val = yplot(midpoint)
yplot = yplot-val
slitmap(*,*,1)=slitmap(*,*,1)-val
xval = xplot(midpoint)
xplot = xplot-xval
badvel = where(yploterr gt 200.) 
if badvel(0) ne -1 then remove, badvel, xplot,yplot,yploterr
print,'A velocity shift of ', val,' km/s has been applied to the map.'

;p_vel = errorplot(xplot,yplot,yploterr,xtitle='distance(kpc)',ytitle='velocity (km/s)',symbol='Diamond',linestyle=6,sym_filled=1,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9],name='Velocity')
;p_Vel.SYM_COLOR ="red"
;p_vel.ERRORBAR_COLOR="red"

writefits,'/scr2/nichal/workspace/slitmaps/kineticmodel/'+name+'_slitmap.fits',slitmap
;make another velfield fits file with corrected velocity (subtract the shift)
newfilename = strsplit(velfile,'.',/extract)
newfilename = path+newfilename(0)+'_kineticmodel.fits'
writefits,newfilename,velocityfield-val

;Overplot with velocitydispersion
makerotmap,name,bestguess_angle,xref,yref,xmid,ymid,x,y,n_slits,dispmap,dispmaperr,something,someerr,slit_size,xplotdisp,yplotdisp,yploterrdisp,1,slit_indices

xplotdisp = xplotdisp-xval ;shift the midpoint to be at 0 kpc
badx=[]
for cc=0,n_elements(xplotdisp)-1 do begin
   test = where(xplot eq xplotdisp(cc)) 
   if test(0) eq -1 then badx=[badx,cc]
   ;remove the points where there is velocity point.
endfor
;stop
if n_elements(badx) ne 0 then remove,badx,xplotdisp,yplotdisp,yploterrdisp

;p_disp = errorplot(xplotdisp,yplotdisp,yploterrdisp,/current,overplot=1,symbol='square',linestyle=6,sym_filled=1,thick=2,name='Velocity dispersion')
;p_disp.SYM_COLOR ="dark green"
;p_disp.ERRORBAR_COLOR="dark green"
;;leg = LEGEND(TARGET=[p_vel,p_disp],position=[0.65,0.27],/normal)
;stop
;p_disp.save,'/scr2/nichal/workspace/output/rotationcurve/kineticmodel/'+name+'_rotcurve.png',BORDER=10, RESOLUTION=100
;p_disp.close

set_plot,'ps'
loadct,2
psname='/scr2/nichal/workspace/output/rotationcurve/kineticmodel/'+name+'_rotcurve.eps'
device, filename = psname,xsize = 10,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,xplot,yplot,yploterr,psym=1,errcolor=70,xmargin=[10,7],xtitle='distance(kpc)',ytitle='velocity (km/s)',type=0,font=0,/nodata,yrange=[min([yplot,yplotdisp])-1.,max([yplot,yplotdisp+yploterrdisp])+1.]
cgplot,xplot,yplot,color=70,psym=14,/overplot
oploterror,xplotdisp,yplotdisp,yploterrdisp,errcolor=10,psym=3
cgplot,xplotdisp,yplotdisp,color=10,psym=15,/overplot,symsize=0.75
device,/close

print, 'Max rotation velosity  is',maxvel,'+/-',maxvelerr,'with angle',bestguess_angle,'degree'
print,'velocity dispersion is', mean(yplotdisp),stddev(yplotdisp)


;Metal gradient

;1) plot the N2 index along major axis and measure the gradient with PP04

;Pick only pixels in the slit along major axis


;if name eq 'CSWA20' then best_angle(0)= best_angle(0)-90.

;stop
metalgradientnew_mar15_fix,bestguess_Angle,xref,yref,xmid,ymid,x,y,n_slits,midindex,N2map,N2maperr,distance,N2_distance,N2_distance_err,radius,N2_radius,N2_radius_err,slit_size,name,fitparam,fitparamerr
distance = distance-midindex*slit_size

gradient=string(fitparam(0),format='(F6.3)')
gradienterr=string(fitparamerr(0),format='(F6.3)')
central_metallicity = string(fitparam(1),format='(F6.3)')
central_metallicityerr=string(fitparamerr(1),format='(F6.3)')

pmetal = errorplot(distance,N2_distance,N2_distance_err,xtitle='distance (kpc)',ytitle='[NII]/Ha',symbol='diamond',sym_color='blue',sym_filled=1,errorbar_color='blue',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9])
pmetal.save,'/scr2/nichal/workspace/output/rotationcurve/kineticmodel/'+name+'_metalN2_PP04.png',BORDER=10,RESOLUTION=100
stop
pmetal.close

pmetal = errorplot(radius,N2_radius,N2_radius_err,xtitle='radius from center(kpc)',ytitle='[NII]/Ha',symbol='diamond',sym_color='blue',sym_filled=1,errorbar_color='blue',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9])
fineradius = findgen(100)/100.*max(radius)
model_y = 10.^((fitparam(0)*fineradius+fitparam(1)-8.90)/0.57)
pmetal = plot(fineradius,model_y,symbol='none',color='green',/current,overplot=1)
;read,xtext,prompt='x location for the legend'
;read,ytext,prompt='y location for the legend'
;t1 = TEXT(xtext,ytext,'gradient=' +gradient+'$\pm$'+gradienterr+' dex/kpc',/DATA, FONT_SIZE=14, FONT_NAME='Helvetica')
;t2 = TEXT(xtext,ytext-0.01,'$12+log(O/H)_{central}$=' +central_metallicity+'$\pm$'+central_metallicityerr+' dex',/DATA, FONT_SIZE=14, FONT_NAME='Helvetica')

pmetal.save,'/scr2/nichal/workspace/output/rotationcurve/kineticmodel/'+name+'_metalgradient.png',BORDER=10,RESOLUTION=100
stop
pmetal.close


;Metal gradient from Maiolino08
;if file_test(m08mapfile) eq 1 then begin

;   M08map = readfits(path+M08mapfile,header)
   ;M08errmap = readfits(path+M08err_mapfile)

;   makerotmap,best_angle(0),xref,yref,xmid,ymid,x,y,n_slits,M08map,something,someerr,slit_size,xplotdisp,yplotdisp,yploterrdisp,1


;;pro makerotmap,angle,xref,yref,xmid,ymid,x,y,n_slits,velocityfield,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,veldisptag

  ; pmetal = errorplot(radius_metal,metal,metal_err,xtitle='radius from center(kpc)',ytitle='log(O/H)+12',symbol='diamond',title=name,sym_color='blue',sym_filled=1,errorbar_color='blue',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9])
   ;pmetal = plot(!x.crange,!x.crange*linparam(1)+linparam(0),symbol='none',color='green',/current,overplot=1)
   ;pmetal.save,'/scr2/nichal/workspace/output/rotationcurve/'+name+'_metalgradnewM08.png'
   ;pmetal.close
   ;stop
;endif

end
