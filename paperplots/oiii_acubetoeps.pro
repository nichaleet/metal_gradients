pro oiii_acubetoeps,acubefile,xc,yc,width,type,scale,name
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

nameout_prime= '/scr2/nichal/workspace/paperplots/'+name+'/'+name
nameout_velmaps='/scr2/nichal/workspace/output/chopped_sourceplane/'+name+'/'+name

;OIII
arraytoeps,acube_chopped(*,*,1),scale,nameout_prime+'_sourceoiii.eps',xtitle,ytitle,'flux(erg/s/cm2)',format='(F6.3)'
writefits,nameout_velmaps+'_sourceoiii.fits',acube_chopped(*,*,1)

;velocity field
arraytoeps,acube_chopped(*,*,0),scale,nameout_prime+'_sourcekinematic_oiii.eps',xtitle,ytitle,'velocity(km/s)',format='(I0)'
writefits,nameout_velmaps+'_sourcekinematic_oiii.fits',acube_chopped(*,*,0)
writefits,nameout_velmaps+'_sourcekinematic_oiii_err.fits',acube_chopped(*,*,8)

;velocity dispersion
;fix the velocity dispersion (subtract the intrinsic dispersion)
imgnew = acube_chopped(*,*,2)
vel = imgnew(where(finite(imgnew)))
velsq = vel^2-50.^2
velsq(where(velsq lt 0.)) = 0.
velnew = sqrt(velsq)
imgnew(where(finite(imgnew))) = velnew

arraytoeps,imgnew,scale,nameout_prime+'_sourceveldisp_oiii.eps',xtitle,ytitle,'velocity dispersion(km/s)',format='(I0)'
writefits,nameout_velmaps+'_sourceveldisp_oiii.fits',imgnew
writefits,nameout_velmaps+'_sourceveldisp_oiii_err.fits',acube_chopped(*,*,10)

;psf
cubesize = size(acube_chopped)
if cubesize(3) eq 13 then writefits,nameout_velmaps+'_sourceoiii_psf.fits',acube_chopped(*,*,12)


;Hbeta
arraytoeps,acube_chopped(*,*,6),scale,nameout_prime+'_sourcehb.eps',xtitle,ytitle,'Hb flux(erg/s/cm2)',format='(F6.3)'
writefits,nameout_velmaps+'_sourcehb.fits',acube_chopped(*,*,6)

;OIII_err
writefits,nameout_velmaps+'_sourceoiii_err.fits',acube_chopped(*,*,9)
;Hbeta_err
writefits,nameout_velmaps+'_sourcehb_err.fits',acube_chopped(*,*,11)

stop
end
