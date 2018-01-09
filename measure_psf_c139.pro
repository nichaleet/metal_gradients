pro xfunc,xpix,A,f,pder
aa=a(0)
b=a(1)
c=a(2)
d=a(3)
e=a(4)
g=a(5)
h=a(6)
f = aa+b*exp(-(xpix-c)^2/(2.*d^2))+e*exp(-(xpix-g)^2/(2.*h^2))
pder=fltarr(n_elements(xpix),n_elements(A))
end

pro measure_psf_c139

restore,'c139_psf.sav'
xflux(where(xflux lt 0.)) = 0.
xflux[0:10]=0.
xflux[23:60]=0.
;xflux yflux
xpix = findgen(n_elements(xflux))
ypix = findgen(n_elements(yflux))

a =[0,0.2,17.5,2.,0.07,22.,0.1]
weight=fltarr(n_elements(xflux))+1.
xfit = curvefit(xpix,xflux,weights,a,sigma,fita=[1,1,1,1,1,1,1],function_name='xfunc',/noderivative)
plot,xpix,xflux
oplot,xfit
stop
end
