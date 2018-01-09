pro c31ana
z=1.5
setplot,14
img  =readfits('c31.fits')

gal1 = img(*,104:117)
gal2 = img(*,15:28)
bg = img(*,29:39)
galarr = [[[gal1]],[[gal2]]]
;gal1
for j=0,1 do begin
   gal=galarr[*,*,j]
   size=size(gal)
   wl=findgen(size[1])
   yerr = fltarr(size[1])
   for i=0, size[1]-1 do yerr(i)=stdev(bg(i,*))

   wl = ((findgen(size[1])-30.)*1.62896+(6563*(z+1)))/10000.
   total_spec = total(gal,2)/size[2]
   ploterror,findgen(size[1]),total_spec,yerr
   ;stop
   
   plot,wl,total_spec,psym=10,xtitle='um',title='CSWA31 Im'+strtrim(j+1,2)
;fit Halpha
   y    = total_spec
   yfit = y(25:39)
   wlfit= wl(25:39)
   measure_error=yerr(25:39)
   fit = gaussfit(wlfit,yfit,A,nterms=4,sigma=sigma,measure_errors=measure_error)

   ;A=fltarr(4)-1.
   ;while A(0) lt 0. do begin
   ;   remove,[n_elements(yfit)-1],wlfit,yfit
   ;   fit = gaussfit(wlfit,yfit,A,nterms=4,sigma=sigma)
   ;endwhile    
   yest = A(0)*exp(-0.5*((wl-A(1))^2/A(2)^2))+A(3)
   oplot,wl,yest,color=50
   area_Ha= A(0)*A(2)*sqrt(2.*!pi)
   area_Ha_err =area_Ha*sqrt((sigma(0)/a(0))^2+(sigma(2)/a(2))^2)

;fit second line
   nii_peak = A(1)+(z+1.)*20.6/10000.
   print,'NII is supposed to be at', nii_peak
   wlmin= nii_peak-4.3*A(2)
   wlmax= nii_peak+4.3*A(2)
   good = where(wl gt wlmin and wl lt wlmax)
   wlNII = wl(Good)
   NII = y(good)
   weightNII=fltarr(n_elements(NII))+1.
   fit2ndline_fixbg,A(2),nii_peak,A(3),wlNII,NII,weightNII,fit,area,sigma_area
   oplot,wlNII,fit,color=20

   NII_Ha = area/area_Ha
   NII_Ha_err = NII_Ha*sqrt((area_Ha_err/area_Ha)^2+(sigma_area/area)^2)
   print,'NII/Ha= ', NII_Ha,NII_Ha_err
   xyouts,400,100,'[NII]/Ha='+strtrim(NII_Ha,2)+'+-'+strtrim(NII_Ha_err,2),/device
   stop
   endfor
end
