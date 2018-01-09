function calymodel_2d_conv,x1,x2,param,kernel,ind,bad,count,name

i  = param(0)    ;inclination
theta = param(1) ;position angle
xc = param(2)    ;center coordinate 
yc = param(3)    ;center coordinate
Rt = param(4)    ;scale radius
vc = param(5)    ;asymptotic velocity
v0 = param(6)    ;systemic velocity


;need to convolve with the 2d gaussian of the psf in source plane image
;make array of the predicted data
sizex=max(ind(0,*))+1
sizey=max(ind(1,*))+1
xarr= rebin(findgen(sizex),sizex,sizey)
yarr= rebin(transpose(findgen(sizey)),sizex,sizey)

x = xarr-xc 
y = yarr-yc
;x and y is now distance to the center in x and y direction in the projected plane

;see Log book page 93
theta_rad = theta/180.*!pi
irad      = i/180.*!pi
znorth    = (x*sin(theta_rad)+y*cos(theta_rad))/sin(irad)
xnorth    = x*cos(theta_rad)-y*sin(theta_rad)   
R         = sqrt(znorth^2+xnorth^2)
VR        = v0+2./!pi*vc*atan(R/Rt)
Vz        = cos(irad)*cos(atan(znorth,xnorth))*VR


;convolve zarr with the 2d gaussian

Vz_conv = convol(Vz,kernel)
im = [[[Vz]],[[Vz_conv]]]
mkhdr,header,im
history = ['i,theta,xc,yc,vc,v0,rt=',string([i,theta,xc,yc,vc,v0,rt])]
sxaddhist,history,header,/comment
ycal=vz_conv
ycal(bad) = 1./0.
;stop
 if (count mod 100000.) eq 0. then writefits,name+'_velmodel.fits',[[[Vz]],[[Vz_conv]],[[ycal]]],header


remove,bad,ycal

return,ycal

end

