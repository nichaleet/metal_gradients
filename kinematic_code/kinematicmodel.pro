pro kinematicmodel,velfile,velerrfile,hadetectionfile,kernel,i0,theta0,xcmin,xcmax,ycmin,ycmax,name,thetamin=thetamin,thetamax=thetamax,imin=imin,imax=imax,vcmax=vcmax,vmeanmin=vmeanmin,vmeanmax=vmeanmax,rtmax=rtmax,pixscale=pixscale
;This program will find the best parameters (inclination i, position angle, disc center coordinate, scale radius, asymptotic velocity, and systemic velocity) that fit the velocity map. (seef Jones, 2010a (the MNRAS one)). Use MCMC method.

;normalize the kernel
if n_elements(kernel) eq 0. then kernel = [1.]
kernel = kernel/total(kernel)

;read the velocity maps
velmap = readfits(velfile,hdr)
velerr = readfits(velerrfile)
velerr(where(velerr eq 0.)) = 5. 

;check the Ha detection
detection = readfits(hadetectionfile)
bad = where(detection le 10.)
velmap(bad) = 1./0.
velerr(bad) = 1./0.

;Keep the original maps
velmap_old=velmap
velerr_old=velerr

if name eq '0744' then begin
   velmap[44,68] = 1./0.
   velerr[44,68] = 1./0.
   velmap[34,71] = 1./0.
   velerr[34,71] = 1./0.
   velmap[57,68] = 1./0.
   velerr[57,68] = 1./0.
   velmap[94,69] = 1./0.
   velerr[94,69] = 1./0.
   velmap[96,63] = 1./0.
   velerr[96,63] = 1./0.
endif

if name eq 'cswa11' then begin
   bad = where(velmap lt -100.)
   velmap(bad)=1./0.
   velerr(bad)=1./0.
   velmap[0:10,0:9] = 1./0.
   velerr[0:10,0:9] = 1./0.
   velmap[*,0:3]=1./0.
   velerr[*,0:3]=1./0.
endif

if name eq 'cswa19' then begin
   velmap[0:31,*]=1./0.
   velerr[0:31,*]=1./0.
   velmap[*,0:10]=1./0.
   velerr[*,0:10]=1./0.
endif

if strmatch(name,'cswa28*') then begin
   bad = where(velmap lt -120. or velmap gt 50.)
   velmap(bad) =1./0.
   velerr(bad) = 1./0.
endif

if name eq 'cswa31' then begin
   velmap[0:5,13:16]=1./0.
   velerr[0:5,13:16]=1./0.
endif

if name eq 'cswa128' then begin
   velmap[0:10,0:6]=1./0.
   velerr[0:10,0:6]=1./0.
   velmap[22:26,26]=1./0.
   velerr[22:26,26]=1./0.
   bad = where(velmap gt 250.)
   velmap(bad)=1./0.
   velerr(bad)=1./0.
endif

if name eq 'cswa139' then begin
   ;make the bottom blob into a merger
   yindex=rebin(transpose(findgen(36)),36,36)
   bad = where(yindex le 15 and velmap gt 30.)
   ;bad = [bad, where(yindex gt 30.)]
   velmap(bad) = 1./0.
   velerr(bad) =1./0.
endif

if name eq 'cswa159' then begin
   hadetection=readfits('/scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps/prettymaps/crop/sourceha_detection_pretty.fits')
   bad=where(hadetection lt 5.5)
   velmap(bad)=1./0.
   velmap_old=velmap
   velerr_old=velerr
   
endif

if name eq 'cswa165' then begin
   bad = where(velmap gt 80.)
   velmap(bad) = 1./0.
   velerr(bad) = 1./0.
   velmap[24,*]=1./0.
   velerr[24,*]=1./0.
endif

if name eq 'abell773' then begin 
;the psf is very elongated. have to expand the size of the velmap
   velmap = fltarr(70,70)+1./0.
   velerr = velmap
   mapsize = size(velmap_old)
   velmap[10:mapsize(1)+9,10:mapsize(2)+9] = velmap_old
   velerr[10:mapsize(1)+9,10:mapsize(2)+9] = velerr_old
   xcmin=xcmin+10
   ycmin=ycmin+10
   xcmax=xcmax+10
   ycmax=ycmax+10
   velmap_old=velmap
   velerr_old=velerr
   newsize = size(velmap,/dimension)
   xind = rebin(findgen(newsize(0)),newsize)
   bad = where(xind lt 22. and velmap lt -60.)
   velmap(bad) = 1./0.
   velerr(bad) = 1./0.
endif

if name eq 'cswa19_kmos' then begin
   velmap[5,8] = 1./0.
   velerr[5,8]=1./0.
   velmap[7,8] = 1./0.
   velerr[7,8]=1./0.
endif

writefits,name+'_Ha_velmap.fits',[[[velmap_old]],[[velerr_old]],[[velmap]],[[velerr]]],hdr


;The two maps must have the same size.
if n_Elements(velmap) ne n_elements(velerr) then stop 

ind = array_indices(velmap,indgen(n_elements(velmap)))
x = ind(0,*)
y = ind(1,*)

nparam = 7
nstep = 200000L

paraname = ['inclination','position angle','center x','center y','scale radius','asymptotic velocity','systemic velocity']
stepsizes = [0.5,1.,1.,1.,.1,3.,3.]*0.5

fixparam= [1,1,1,1,1,1,1] ;0 for params you want to fix
if xcmin eq xcmax then fixparam(2) =0.
if ycmin eq ycmax then fixparam(3) =0.
stepsizes = stepsizes*fixparam


;guess initial parameters
xc0 = mean([xcmin,xcmax])
yc0 = mean([ycmin,ycmax])
Rt0 = 1.  ;3 pixels
vc0 = 100.
if name eq '0744' then begin
   vc0=600.
   rt0=10.
endif
vmean0 = 0.;velmap(xc0,yc0) ;v at center is the systemic velocity
p0 = [i0,theta0,xc0,yc0,rt0,vc0,vmean0]

;Range of the parameters
if n_elements(imin) eq 0. then imin = 0.
if n_elements(thetamin) eq 0. then thetamin = 0.
if n_elements(vmeanmin) eq 0. then vmeanmin = -10.
rtmin = 0.
vcmin = 10.
rangemin = [imin,thetamin,xcmin,ycmin,rtmin,vcmin,vmeanmin]

if n_elements(imax) eq 0. then imax = 90.
if n_elements(thetamax) eq 0. then thetamax=360.
if n_Elements(vcmax) eq 0. then vcmax=300.
if n_elements(vmeanmax) eq 0. then vmeanmax = 10.
if n_Elements(rtmax) eq 0. then rtmax = 5.

rangemax= [imax,thetamax,xcmax,ycmax,rtmax,vcmax,vmeanmax]

;clean data (get rid of the Nan pixels)
bad = where(finite(velmap) eq 0)
remove,bad,x,y,velmap,velerr

mcmc_2d,x,y,velmap,velerr,nparam,nstep,p0,stepsizes,rangemin,rangemax,paraname,kernel,ind,bad,name,pixscale,fixparam,returnvalues
          ;Return values is an array of size 3*nparam contain 
          ;lower limit,best estimate, and upper limit values for each parameter.

pbest = returnvalues(1,*)
velmodel = calymodel_2d_conv(x,y,pbest,kernel,ind,bad,0,name)
ngoodparam  = n_elements(where(fixparam ne 0.))
chisquarebest = calchisq(velmodel,velmap,velerr,ngoodparam)
residuals = velmodel-velmap
plothist,residuals,title='velocity residuals'
rms_residuals = sqrt(total((velmodel-velmap)^2)/n_elements(velmodel))
print,'RMS of velocity residual =', rms_residuals
print,'best reduced  chisquare =',chisquarebest

if strmatch(name,'cswa28*') then begin
   openw,1,'cswa28result.dat',/append
   printf,1,'RMS of velocity residual =', rms_residuals
   printf,1,'best reduced  chisquare =',chisquarebest
   printf,1,'RT (kpc)',returnvalues(*,4)*0.424
   close,1
endif

end
