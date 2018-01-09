pro makerotmap,angle,xref,yref,x,y,n_slits,velocityfield,rotationvelocity,slit_size
  angle_perpend = angle-90.
  slope = tan(angle_perpend/180.*!pi)

; y = slope*x+intercept
; intercept = y-slope*x
  intercept = yref-slope*xref

  velocity = fltarr(n_slits)
  velocity_err = fltarr(n_slits)
  for i=0,n_slits-1 do begin
     pix_in_slit = where(y le slope*x+intercept(i+1) and y ge slope*x+intercept(i))
     goodvel = velocityfield(pix_in_slit)
     if pix_in_slit(0) ne -1 then begin 
        print,'There are', n_Elements(where(finite(goodvel) eq 1)),'pixels with velocity in this slit'
        velocity(i) = mean(goodvel(where(finite(goodvel) eq 1)))
        velocity_err(i) = stddev(goodvel(where(finite(goodvel) eq 1)))
     endif else print,'There is no pixel in this slit.'
  endfor
  ploterror,findgen(n_slits)*slit_size,velocity,velocity_err,xtitle='radius(kpc)',ytitle='velocity (km/s)',psym=1
  rotationvelocity = max(velocity(where(finite(velocity) eq 1)))-min(velocity)
  ;stop
end

pro rotationcurve,Hafile,xc,yc,width,bestguess_angle,pixelscale,z,metalmap
setplot,14
;example for cswa128
;hafile='/scr2/nichal/workspace/output/cswa128_Ha_Kn2_mosaic_sky_130hr_acube_sourceplane_interp.fits'
;x=245  ; center of the image (approx)
;y=237
;width = 24.  ;width in pixel.should be an even number.
;bestguess = 135. ;angle of major axis in degree.

halfwidth = width/2.
slitwidth =  2.  ;slit width in pixel
n_slits = fix(width/slitwidth)

;Calculate physical size per pixel
theta_1kpc = zang(1.,z) ; in arcsecond
pix_size  = pixelscale/theta_1kpc ;in kpc per pixel

Hacube = readfits(hafile,header)

if hafile eq '/scr2/nichal/workspace/output/cswa20_Ha_Hn2_mosaic_scaledsky_100hr_acube_sourceplane_interp.fits' then hacube(1,28,*)=1./0.

velocityfield = hacube(xc-halfwidth:xc+halfwidth,yc-halfwidth:yc+halfwidth,0) ; size is (width+1,width+1)

;Change the center of the field from x,y original to the coordinate of the cropped field where the bottom left is (0,0) which is the x,y position of the parameter velocityfield

xmid = halfwidth
ymid = halfwidth

xref = fltarr(n_slits+1)+xmid  ;xref, yref are the points that each slits are defined
yref = findgen(n_slits+1)*slitwidth

ind = array_indices(velocityfield,findgen(n_elements(velocityfield)-1))
x = ind(0,*)
y = ind(1,*)


; Do for loop for n degrees around the best guess angle
nang = 30
maxvel = fltarr(nang+1)
angdiff = findgen(nang+1)-nang/2.
angle = bestguess_angle+angdiff

for angind = 0,nang do begin
anglenow = angle(angind)
slit_size = slitwidth*pix_size*cos((anglenow-90.)*!pi/180.)
makerotmap,anglenow,xref,yref,x,y,n_slits,velocityfield,rotationvelocity,slit_size
maxvel(angind) = rotationvelocity
;wait,1.
endfor

plot, angle, maxvel
best_angle = angle(where(maxvel eq max(maxvel)))
print, 'Max rotation velosity  is',max(maxvel),'with angle',best_angle,'degree'
;stop

makerotmap,best_angle(0),xref,yref,x,y,n_slits,velocityfield,rotationvelocity,slit_size

stop

;Metal gradient
metalcube = readfits(metalmap,header)
metalMaiolino  = metalcube(xc-halfwidth:xc+halfwidth,yc-halfwidth:yc+halfwidth,2)
metalgradient,metalmaiolino,xmid,ymid,x,y,n_slits,slit_size/2.,slitwidth/2.
stop
end
