pro makerotmap,angle,xref,yref,x,y,n_slits,velocityfield,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,veldisptag
  
  angle_perpend = angle-90.
  slope = tan(angle_perpend/180.*!pi)

; y = slope*x+intercept
; intercept = y-slope*x
  intercept = yref-slope(0)*xref

  velocity = fltarr(n_slits)
  velocity_err = fltarr(n_slits)
  for i=0,n_slits-1 do begin
     pix_in_slit = where(y le slope*x+intercept(i+1) and y ge slope*x+intercept(i))
     goodvel = velocityfield(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(goodvel) eq 1)),'pixels with velocity in this slit'
        velocity(i) = mean(goodvel(where(finite(goodvel) eq 1)))
        velocity_err(i) = stddev(goodvel(where(finite(goodvel) eq 1)))
;        if veldisptag eq 0 then begin
;           velocity(i) = mean(goodvel(where(finite(goodvel) eq 1)))
;           velocity_err(i) = stddev(goodvel(where(finite(goodvel) eq 1)))
;        endif else begin
;           goodvel = goodvel(where(finite(goodvel)eq 1))
;           goodvel = sqrt(goodvel^2-50.^2) ;Minus the instrument intrinsic dispersion
;           velocity(i) = mean(goodvel)
;           velocity_err(i) = stddev(goodvel)
;        endelse

     endif else print,'There is no pixel in this slit.'
  endfor
  xplot = findgen(n_slits)*slit_size

  badpix=where(finite(velocity_err) eq 0)
  if badpix(0) ne -1 then remove, xplot,velocity,velocity_err

  yplot=velocity
  yploterr=velocity_err
  
  ploterror,xplot,velocity,velocity_err,xtitle='radius(kpc)',ytitle='velocity (km/s)',psym=1
  currentmaxvel = max(velocity(where(finite(velocity) eq 1)))

  rotationvelocity = currentmaxvel-min(velocity)
  velerr1 = velocity_err(where(velocity eq currentmaxvel))
  velerr2 = velocity_err(where(velocity eq min(velocity)))
  rotationvelocityerr = sqrt(velerr1^2+velerr2^2)
  rotationvelocityerr = rotationvelocityerr(0)



 ; stop
end

pro rotationcurve,velfile,xc,yc,bestguess_angle,pixelscale,metalmap,name

;example for cswa128
;rotationcurve,'/scr2/nichal/workspace/output/cropvel/cswa128.fits',14,14,135.,0.04,2.22,'/scr2/nichal/workspace/output/cropmetal/cswa128.fits'
;pixelscale is in kpc per pixel

velocityfield= readfits(velfile)
if name eq 'CSWA20' then velocityfield(where(velocityfield ge 120.)) = 1./0.
;stop
sizeim = size(velocityfield)
width = sizeim(1)
slitwidth =  2.  ;slit width in pixel
n_slits = fix(width/slitwidth)

;Calculate physical size per pixel
;theta_1kpc = zang(1.,z) ; in arcsecond
;pix_size  = pixelscale/theta_1kpc ;in kpc per pixel
pix_size = pixelscale ; in kpc per pixel


xmid = xc
ymid = yc

xref = fltarr(n_slits+1)+xmid  ;xref, yref are the points that each slits are defined
yref = findgen(n_slits+1)*slitwidth

ind = array_indices(velocityfield,findgen(n_elements(velocityfield)-1))
x = ind(0,*)
y = ind(1,*)


; Do for loop for nang degrees around the best guess angle
nang = 30
maxvel = fltarr(nang+1)
maxvelerr = fltarr(nang+1)
angdiff = findgen(nang+1)-nang/2.
angle = bestguess_angle+angdiff

for angind = 0,nang do begin
anglenow = angle(angind)
slit_size = slitwidth*pix_size*cos((anglenow-90.)*!pi/180.)
set_plot,'x'

makerotmap,anglenow,xref,yref,x,y,n_slits,velocityfield,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,0
maxvel(angind) = rotationvelocity
maxvelerr(angind) = rotationvelocityerr
;wait,1.
endfor

plot, angle, maxvel
good = where(finite(maxvel) eq 1)
angle = angle(good)
maxvel = maxvel(good)
maxvelerr = maxvelerr(good)
best_angle = angle(where(maxvel eq max(maxvel)))
print, 'Max rotation velosity  is',max(maxvel),'+/-',maxvelerr(where(maxvel eq max(maxvel))),'with angle',best_angle,'degree'


;Plot the best rotation curve to file output
setplot,39
makerotmap,best_angle(0),xref,yref,x,y,n_slits,velocityfield,rotationvelocity,rotationvelocityerr,slit_size,xplot,yplot,yploterr,0

file_mkdir,'/scr2/nichal/workspace/output/rotationcurve'
pos1 = strpos(velfile,'/',/reverse_Search)
pos2 = strpos(velfile,'.fits')
fileout = '/scr2/nichal/workspace/output/rotationcurve/'+strmid(velfile,pos1,pos2-pos1)+'_rotationcurve.eps'
loadct,39

;shift the median velocity to 0
medpoint = n_elements(yplot)/2
val = yplot(medpoint)
yplot = yplot-val
p_vel = errorplot(xplot,yplot,yploterr,xtitle='radius(kpc)',ytitle='velocity (km/s)',title=name,symbol='Diamond',linestyle=6,sym_filled=1,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9])
p_Vel.SYM_COLOR ="red"
p_vel.ERRORBAR_COLOR="red"

;Overplot with velocitydispersion
dispfile = strmid(velfile,0,pos2)+'_disp.fits'
dispmap = readfits(dispfile)
makerotmap,best_angle(0),xref,yref,x,y,n_slits,dispmap,something,someerr,slit_size,xplotdisp,yplotdisp,yploterrdisp,1

;if name eq 'cswa128' then remove,[0],xplotdisp,yplotdisp,yploterrdisp
;plot,xplotdisp,yplotdisp,color=50
p_disp = errorplot(xplotdisp,yplotdisp,yploterrdisp,/current,overplot=1,symbol='square',linestyle=6,sym_filled=1,thick=2)
p_disp.SYM_COLOR ="dark green"
p_disp.ERRORBAR_COLOR="dark green"

p_disp.save,'/scr2/nichal/workspace/output/rotationcurve/'+name+'_rotcurve.png'
p_disp.close

print,'velocity dispersion is', mean(yplotdisp),stddev(yplotdisp)

;Metal gradient
metalcube = readfits(metalmap,header)
metalgradient,metalcube,xmid,ymid,x,y,n_slits,slit_size/2.,slitwidth/2.,radius_metal,metal,metal_err,linparam
fileout = '/scr2/nichal/workspace/output/rotationcurve/'+strmid(velfile,pos1,pos2-pos1)+'_metalN2.png'
pmetal = errorplot(radius_metal,metal,metal_err,xtitle='radius from center(kpc)',ytitle='log(O/H)+12',symbol='diamond',title=name,sym_color='blue',sym_filled=1,errorbar_color='blue',linestyle=6,thick=2,font_size=20,position=[0.2,0.15,0.95,0.9])
pmetal = plot(!x.crange,!x.crange*linparam(1)+linparam(0),symbol='none',color='green',/current,overplot=1)
pmetal.save,'/scr2/nichal/workspace/output/rotationcurve/'+name+'_metalgrad.png'
pmetal.close


;Metal gradient from Maiolino08
metalmap = strmid(metalmap,0,strpos(metalmap,'.fits'))+'m08.fits'

print, 'Maiolino08 gradient'
metalcube = readfits(metalmap,header)
metalgradient,metalcube,xmid,ymid,x,y,n_slits,slit_size/2.,slitwidth/2.,radius_metal,metal,metal_err,linparam
;stop
end
