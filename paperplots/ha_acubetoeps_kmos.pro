pro ha_acubetoeps,acubefile,xc,yc,width,type,scale,name
;type = 'kpc' then the output eps will have dimension in kpc
;type = 'arcsec' then the output eps will have dimension in arcsec
;scale = size per pixel according to type
;name = name of this galaxy
;output will be in /workspace/paperplots

acube=readfits(acubefile,hdr)
hwidth = fix(width/2)

lowx=xc-hwidth
upx =xc+hwidth
lowy=yc-hwidth
upy =yc+hwidth
sz = size(acube)

if lowx lt 0. then lowx=0
if lowy lt 0. then lowy=0
if upx ge sz(1) then upx=sz(1)-1
if upy ge sz(2) then upy=sz(2)-1

if hwidth ne 0. then acube_chopped = acube(lowx:upx,lowy:upy,*) else acube_chopped=acube


ytitle=name
xtitle = type

file_mkdir,'/scr2/nichal/workspace/paperplots/'+name
file_mkdir,'/scr2/nichal/workspace/output/chopped_sourceplane/'+name

nameout_prime= '/scr2/nichal/workspace/paperplots/'+name+'/'+name
nameout_velmaps='/scr2/nichal/workspace/output/chopped_sourceplane/'+name+'/'+name

;velocity field
arraytoeps,acube_chopped(*,*,0),scale,nameout_prime+'_sourcekinematic.eps',xtitle,ytitle,'velocity(km/s)',format='(I0)'
writefits,nameout_velmaps+'_sourcekinematic.fits',acube_chopped(*,*,0)
writefits,nameout_velmaps+'_sourcekinematic_err.fits',acube_chopped(*,*,8)

;Ha_Detection
arraytoeps,acube_chopped(*,*,4),scale,nameout_prime+'_sourceha_detection.eps',xtitle,ytitle,'Ha flux(erg/s/cm2)',format='(F6.3)'
writefits,nameout_velmaps+'_sourceha_detection.fits',acube_chopped(*,*,4)


;Halpha
arraytoeps,acube_chopped(*,*,1),scale,nameout_prime+'_sourceha.eps',xtitle,ytitle,'Ha flux(erg/s/cm2)',format='(F6.3)'
writefits,nameout_velmaps+'_sourceha.fits',acube_chopped(*,*,1)


;velocity dispersion
;fix the velocity dispersion (subtract the intrinsic dispersion)
imgnew = acube_chopped(*,*,2)
vel = imgnew(where(finite(imgnew)))
velsq = vel^2-50.^2
velsq(where(velsq lt 0.)) = 0.
velnew = sqrt(velsq)
imgnew(where(finite(imgnew))) = velnew

arraytoeps,imgnew,scale,nameout_prime+'_sourceveldisp.eps',xtitle,ytitle,'velocity dispersion(km/s)',format='(I0)'
writefits,nameout_velmaps+'_sourceveldisp.fits',imgnew
writefits,nameout_velmaps+'_sourceveldisp_err.fits',acube_chopped(*,*,10)

;psf
cubesize = size(acube_chopped)
if cubesize(3) ge 13 then writefits,nameout_velmaps+'_sourceha_psf.fits',acube_chopped(*,*,12)

;continuum
if cubesize(3) ge 14 then writefits,nameout_velmaps+'_sourcecontinuum.fits',acube_chopped(*,*,13)

;NII
arraytoeps,acube_chopped(*,*,6),scale,nameout_prime+'_sourcenii.eps',xtitle,ytitle,'NII flux(erg/s/cm2)',format='(F6.3)'
writefits,nameout_velmaps+'_sourcenii.fits',acube_chopped(*,*,6)

;stop

end
