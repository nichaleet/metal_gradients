pro metalnew_acubetoeps,acubefile,xc,yc,width,type,scale,name
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

;[NII]/Ha
arraytoeps,acube_chopped(*,*,0),scale,nameout_prime+'_sourcen2index.eps',xtitle,ytitle,'[NII]/Ha',format='(F5.2)'
writefits,nameout_velmaps+'_sourcen2index.fits',acube_chopped(*,*,0)
writefits,nameout_velmaps+'_sourcen2index_err.fits',acube_chopped(*,*,3)

;12+log(o/h) M08
arraytoeps,acube_chopped(*,*,2),scale,nameout_prime+'_sourceBayesianmetal.eps',xtitle,ytitle,'12+log[O/H](M08)',format='(F5.2)'
writefits,nameout_velmaps+'_sourceBayesianmetal.fits',acube_chopped(*,*,2)
writefits,nameout_velmaps+'_sourceBayesianmetal_err.fits',acube_chopped(*,*,5)

;O3N2 index
arraytoeps,acube_chopped(*,*,1),scale,nameout_prime+'_sourceo3n2index.eps',xtitle,ytitle,'O3N2 index',format='(F5.2)'
writefits,nameout_velmaps+'_sourceo3n2index.fits',acube_chopped(*,*,1)
writefits,nameout_velmaps+'_sourceo3n2index_err.fits',acube_chopped(*,*,4)

;Metal type
arraytoeps,acube_chopped(*,*,6),scale,nameout_prime+'_sourcemetaltype.eps',xtitle,ytitle,'Metal type',format='(F5.1)'
writefits,nameout_velmaps+'_sourcemetaltype.fits',acube_chopped(*,*,6)
     

end
