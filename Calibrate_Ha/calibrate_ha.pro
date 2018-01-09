pro calibrate_ha,ttcubefile,name=name,filter=filter,photmag=photmag,z=z,calibfactor=calibfactor
;Input: ttcubefile: cube file of tip tilt star
;       photmag   : photometric magnitude of the tip/tilt star. Assume 2MASS magnitude


;read file. Calculate wavelength
starcube = readfits(ttcubefile,hdr)
wl       = getwl_filter(filter) ; wavelength in microns

;remove the obvious glitches
glitches = where(abs(starcube) gt 10.,ctglitch)
if ctglitch gt 1 then starcube(glitches) = median(starcube)

;collapse the cube in wavelength to see the psf
star_2d = total(starcube,1)*0.
star_2d_sigma = star_2d
imsize = size(star_2d,/dimensions)
for ii=0,imsize(0)-1 do begin
   for jj=0,imsize(1)-1 do begin
      meanclip,starcube[*,ii,jj],meanflux,sigma
      star_2d(ii,jj) = meanflux
      star_2d_sigma(ii,jj) = sigma
   endfor
endfor

;fix the blank pixels
star_2d_sigma(where(star_2d eq 0.)) = mean(star_2d(where(star_2d ne 0.)))
star_2d(where(star_2d eq 0.)) = mean(star_2d(where(star_2d ne 0.)))

if name eq 'cswa15' then begin
   xc=32
   yc=19
   peakval=0.04
   goto, skip
endif
if name eq 'cswa28' then begin
   star_2d = readfits('cswa28_2d_ha.fits')
   xc      = 32
   yc      = 17
   peakval = 20
   goto,skip
endif
if name eq 'cswa31' then begin
   xc = 34
   yc = 27
   peakval = 0.2
   goto,skip
endif
if name eq 'cswa165' then begin
   xc = 47
   yc = 25
   peakval = 0.1
   goto,skip
endif
;Find stars
find,star_2d,x,y,flux,sharp,round,0.04,4.,[-0.5,0.5],[0.2,1.]
real_star_pos = where(flux eq max(flux))
xc=(x(real_star_pos))[0]
yc=(y(real_star_pos))[0]
peakval = star_2d(xc,yc)

skip:
;Fit 2d Gaussian
int_param = [0.,peakval,2.,2.,xc,yc,0.]
pi = replicate({fixed:0},7)
pi([4,5]).fixed = 1  ;fix the central locations

psf = mpfit2dpeak(star_2d,param,error=star_2d_sigma,estimates=int_param,parinfo=pi,status=status)
print,'position of star:',xc,yc
print,'FWHM(pixels):',param([2,3])*2.35482
FWTM = round(param([2,3])*4.29193)
if name eq 'cswa28'  then FWTM = round(param([2,3]))
if name eq 'cswa31' then FWTM = round(param([2,3])*1.9)
if name eq 'cswa139' then FWTM = round(param([2,3])*2.5)
if name eq 'cswa165' then FWTM = round(param([2,3])*3.)
if name eq 'abell773' then FWTM[0]=8
print,'FWTM(pixels):',param([2,3])*4.29193,FWTM

;make star region mask
xarr = rebin(findgen(imsize(0)),imsize(0),imsize(1))-xc
yarr = rebin(transpose(findgen(imsize(1))),imsize(0),imsize(1))-yc
   ;Elipse equation: total(distance from the point to each focus) = 2a
a = float(max(FWTM))
b = float(min(FWTM))
c = sqrt(a^2-b^2)
if FWTM(0) gt FWTM(1) then rarr = sqrt((xarr+c)^2+yarr^2)+sqrt((xarr-c)^2+yarr^2) $
else rarr = round(sqrt((yarr+c)^2+xarr^2)+sqrt((yarr-c)^2+xarr^2))

mask = fltarr(imsize(0),imsize(1))
mask(where(rarr le 2.*a)) = 1.

;Plots
zoomfactor = 5.
imsizenew = imsize*zoomfactor
star_2d_zoom = congrid(star_2d,imsizenew(0),imsizenew(1),/interp)
window,0,xsize=imsizenew(0),ysize=imsizenew(1)
tvscl,star_2d_zoom

;Get flux at Ha wavelength in num_res resolution element
      ;The unit in reduced data cube is DN/sec
      ;The spectral resolution is about 3000 -> ~100 km/s
num_Res = 3.
Ha_wl   = (z+1.)*.656461
wl_pix  = sxpar(hdr,'CDELT1')/1000. ;microns
if wl_pix-(wl(1)-wl(0)) gt 1.e-10 then stop,'Stop: Problems with wavelength resolution.'
pix_res = num_Res*round(Ha_wl/3000./wl_pix) 
if pix_res mod 2. ne 1. then pix_res = pix_res+1.
minval = min(abs(wl-ha_wl),pos_ha)
pos_min = pos_ha-0.5*(pix_res-1.)
pos_max = pos_ha+0.5*(pix_res-1.)
print,'pixel of Ha',pos_min,pos_ha,pos_max
      ;collapse in wavelength
star_2d_Haflux_simple = total(starcube[pos_min:pos_max,*,*],1)
star_2d_Haflux_meanclipsum = fltarr(imsize(0),imsize(1))
star_2d_Haflux_meanclipsum_sigma = fltarr(imsize(0),imsize(1))
for ii=0,imsize(0)-1 do begin 
   for jj=0,imsize(1)-1 do begin
      meanclip,starcube[pos_min:pos_max,ii,jj],meanha,sigmaha
      star_2d_Haflux_meanclipsum[ii,jj] = meanha*pix_res
      star_2d_Haflux_meanclipsum_sigma[ii,jj] = sigmaha*pix_res
   endfor
endfor

;write fits
writefits,name+'_star_2d.fits',[[[star_2d]],[[star_2d_Haflux_meanclipsum]],[[mask]],[[rarr]],[[psf]],[[star_2d_sigma]]]

     ;collapse in spatial area
star_reg = where(mask eq 1)
haflux_simple   = total(star_2d_haflux_simple(star_reg))/(pix_res*wl_pix)    ;ADU/second/micron/area
haflux_meanclip = total(star_2d_haflux_meanclipsum(star_reg))/(pix_res*wl_pix) ;ADU/second/micron/area
haflux_meanclip_sigma = sqrt(total(star_2d_Haflux_meanclipsum_sigma^2))
     ;if the two methods differ more than 10% then stop
if abs(haflux_simple-haflux_meanclip)/haflux_meanclip gt 0.1 then stop,'Problems with cosmic rays rejection'
print, haflux_simple,haflux_meanclip

;Calculate Real Flux from stars
    ;what band is the observation?
mainband = strmid(filter,0,1)
case mainband of
  'K': mag_zero_flux = 4.283e-7
  'H': mag_zero_flux = 1.133e-6
  'J': mag_zero_flux = 3.129e-6
  else: print,'No zero flux for this band'
endcase

realflux = mag_zero_flux/(10.^(photmag/2.5)) ;erg/s/cm^2/micron/area

calibfactor = realflux/haflux_meanclip
calibfactor_sigma = calibfactor*haflux_meanclip_sigma/haflux_meanclip
calibfactor = [calibfactor,calibfactor_sigma]
print,name,' conversion_factor = ',calibfactor
stop
end
