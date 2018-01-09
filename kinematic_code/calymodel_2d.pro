function calymodel_2d,x1,x2,param

i  = param(0)    ;inclination
theta = param(1) ;position angle
xc = param(2)    ;center coordinate 
yc = param(3)    ;center coordinate
Rt = param(4)    ;scale radius
vc = param(5)    ;asymptotic velocity
v0 = param(6)    ;systemic velocity


x = abs(x1-xc)
y = abs(x2-yc)

theta_deg = theta/180.*!pi
Rsquare   = (x*cos(theta_deg)-y*sin(theta_deg))^2+((x*sin(theta_deg)+y*cos(theta_deg))/sin(i/180.*!Pi))^2
R         = sqrt(Rsquare)
ycal      = v0+2./!pi*vc*atan(R/Rt)

return,ycal

end

