pro c164ana
z=2.518
setplot,14
img  =readfits('c164.fits')
bg = [[img(*,0:16)],[img(*,72:88)]]

galup = img(*,56:70)
galdown = img(*,17:31)

galarr = [[[galup]],[[galdown]]]

for j=0,1 do begin
   gal=galarr(*,*,j)
   size=size(gal)
   wl=((findgen(size[1])-6.)*2.17+(6562.8*(z+1)))/10000. ;um
   yerr = fltarr(size[1])
   for i=0, size[1]-1 do yerr(i)=stdev(bg(i,*))
   ;yerr(where(yerr ge 100)) = 350.

;align the pixels according to the Ha_wl
   total_spec = total(gal,2)/size[2]
;remove skyline
   total_spec(24:31)=mean(img(0:23,0:21))
   ploterror,findgen(size[1]),total_spec,yerr
   ;stop
   plot,wl,total_spec,psym=10,xtitle='um',title='CSWA164 Im'+strtrim(j+1,2)
;fit Halpha
   y    = total_spec
   yfit = y(5:30)
   wlfit= wl(5:30)
   measure_error = yerr(5:30)
   fit = gaussfit(wlfit,yfit,A,nterms=4,sigma=sigma,measure_errors=measure_error)
   yest = A(0)*exp(-0.5*((wl-A(1))^2/A(2)^2))+A(3)
   oplot,wl,yest,color=50
   area_Ha= A(0)*A(2)*sqrt(2.*!pi)
   area_Ha_err =area_Ha*sqrt((sigma(0)/a(0))^2+(sigma(2)/a(2))^2)
   

;fit second line
   weight = 1./yerr^2
   nii_peak = A(1)+(z+1.)*20.61/10000.
   print,'NII is supposed to be at', nii_peak
   wlmin= nii_peak-6.*A(2)
   wlmax= nii_peak+6.*A(2)
   good = where(wl gt wlmin and wl lt wlmax)
   wlNII = wl(Good)
   NII = y(good)
   weightNII=yerr(good)
   fit2ndline_fixbg,A(2),nii_peak,A(3),wlNII,NII,weightNII,fit,area,sigma_area
   


   oplot,wlNII,fit,color=20

   NII_Ha = area/area_Ha
   NII_Ha_err = NII_Ha*sqrt((area_Ha_err/area_Ha)^2+(sigma_area/area)^2)
   SN=sqrt(abs(total(nii-A(3))))
   print,'SN of NII ~ ',SN
   print,'NII/Ha= ', NII_Ha,NII_Ha_err
   xyouts,200,100,'[NII]/Ha='+strtrim(NII_Ha,2)+'+-'+strtrim(NII_Ha_err,2),/device
   
   stop
endfor

end
