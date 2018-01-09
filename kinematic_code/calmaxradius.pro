function calmaxradius,x1,x2,param,kernel,ind,bad,count,pixscale
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

remove,bad,R

maxradius = pixscale*max(R)

return,maxradius

end
