pro c159reduced_ana
z=2.30
setplot,14
img  =readfits('c159_reduced.fits')
gal_neg = img(*,1003:1010)
gal_pos = img(*,989:996)
bg = [[img(*,980:988)],[img(*,1015:1025)]]

galarr = [[[gal_pos]],[[gal_neg]]]
size=size(gal_pos)
wlfull=((findgen(size[1])-1124.)*2.17195+21651.267)/10000. ;um
yerrfull = fltarr(size[1])
for i=0, size[1]-1 do yerrfull(i)=stdev(bg(i,*))*sqrt(size[2])
;yerr(where(yerr ge 20)) = 50.

wl = wlfull(1100:1170)
yerr = yerrfull(1100:1170)
total_spec_arr=fltarr(n_elements(wl),2)
!p.multi=[0,1,2]
for j=0,1 do begin
   ;j=0 is for positive, j=1 is for negative
   gal=galarr[1100:1170,*,j]
   total_spec_arr(*,j) = total(gal,2)
   ploterror,findgen(n_elements(wl)),total(gal,2),yerr
endfor
stop

!p.multi=[0,1,1]
total_spec = total_spec_arr(*,0)-total_spec_arr(*,1)

plot,wl,total_spec,psym=10,xtitle='um',title='CSWA159'
;fit Halpha
y    = total_spec
yfit = y(20:45)
wlfit= wl(20:45)
measure_error = yerr(20:45)
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
print,'NII/Ha= ', NII_Ha,NII_Ha_err
xyouts,400,100,'[NII]/Ha='+strtrim(NII_Ha,2)+'+-'+strtrim(NII_Ha_err,2),/device
stop
end
