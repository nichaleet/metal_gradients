pro singlet, x, par, f, pder

  x0 = par[0]
  aa = par[1]
  ww = par[2]
  cont = par[3]

  f = (aa/ww/sqrt(2.0*!pi)) * exp(-0.5*(x-x0)^2/ww^2) + cont

  pder = fltarr(n_elements(x),n_elements(par)) ; no value returned.

end

pro fit2ndline,width,center,x,y,weight,fit,area,sigma_area
  a = [center, 2.*width,width,mean(y)]
  fit_a = [0,1,1,1]             ; fit only the intensity and continuum
  fit = curvefit(x,y,weight,a,sigmafit,chisq=chi2, fita=fit_a,function_name="singlet",/noderivative)
  print, 'fit parameter',a
  chi2 = total(weight*(y - fit)^2) ; chi-square of line fit
  chi0 = total(weight*(y-a(3))^2)
  chi2_var = chi0+1. - chi2     ; variance of the line fit
  nsigma = sqrt(chi2_var)       ; significance of line fit, in standard devs
  print, nsigma
  area=a(1)
  sigma_area = area*sqrt((sigmafit(1)/a(1))^2+(sigmafit(2)/a(2))^2)
  width = a(2)
  ;stop
end
