pro bptanalysis,name,OIII,hbeta,NII,halpha,OIII_err,hbeta_err,NII_err,halpha_err

;======================================================
setplot,14
;find detected NII pixels
NII_upper = NII_err+NII
change2upper = where(NII lt 0. and NII_upper gt 0.)
sizeHa = size(halpha)
upper_tags = fltarr(sizeHa(1),sizeHa(2))
upper_tags(where(finite(NII))) = 8.
if change2upper(0) ne -1 then begin
   NII(change2upper) = NII_upper(change2upper)
   upper_tags(change2upper) = 4.
endif
nan = where(NII_upper lt 0.)
upper_tags(nan) = 2.
NII(nan) = 1./0.

;Now tag=8 is detection, tag =4 is upper limit detection, 
                                ;tag=2 is non detection

;find detected Hbeta pixels

Hbeta_upper = Hbeta_err+Hbeta
change2upper = where(Hbeta lt 0. and Hbeta_upper gt 0.)
HbetaUpper_tags = fltarr(sizeHa(1),sizeHa(2))      
HbetaUpper_tags(where(finite(Hbeta))) = 27.

if change2upper(0) ne -1 then begin
   Hbeta(change2upper) = Hbeta_upper(change2upper)
   HbetaUpper_tags(change2upper) = 9.
endif
nan = where(Hbeta_upper lt 0.)
if nan(0) ne -1 then HbetaUpper_tags(nan) = 3.
Hbeta(nan) = 1./0.
;Now tag=27 is detection, tag =9 is upper limit detection, 
                                ;tag=3 is non detection

;NII/Halpha
NII_Ha = alog10(NII/halpha)
NII_Ha_err = sqrt((NII_err/NII)^2+(Halpha_err/Halpha)^2)/alog(10.)
;sigma(logA) = sigma(a)/A/log10
;OIII/Hbeta
OIII_Hb = alog10(OIII/Hbeta)
OIII_Hb_err = sqrt((OIII_err/OIII)^2+(Hbeta_err/Hbeta)^2)/alog(10.)

;make a plot
;1) for all detections
goodnow = where(upper_tags eq 8. and hbetaUpper_tags eq 27. and halpha ge 0.025)
plotsym,8,1,color=200,/fill
ploterror,NII_Ha(goodnow),OIII_Hb(goodnow),NII_Ha_err(goodnow),OIII_Hb_err(goodnow),xtitle='log([NII]6584/Ha)',ytitle='log([OIII]5007/Hb',psym=8,xrange=[-2,0],yrange=[-0.5,1.5],errcolor=200,errthick=1

goodnow = where(upper_tags eq 8. and hbetaUpper_tags eq 27. and halpha le 0.025)
plotsym,0,1,color=30,/fill
oploterror,NII_Ha(goodnow),OIII_Hb(goodnow),NII_Ha_err(goodnow),OIII_Hb_err(goodnow),errcolor=30,errthick=1,psym=8

;2) for NII upper limit
goodnow = where(upper_tags eq 4 and hbetaupper_tags eq 27)
if goodnow(0) ne -1 then begin
   plotsym,6,2,color=60
   oplot,NII_Ha(goodnow),OIII_Hb(goodnow),psym=8
   print, 'There are',n_elements(goodnow),'NII upper limit'
endif
;3) for Hb upper limit
goodnow = where(upper_tags eq 8 and hbetaupper_tags eq 9)
if goodnow(0) ne -1 then begin
   plotsym,1,2,color=150
   oplot,NII_Ha(goodnow),OIII_Hb(goodnow),psym=8
   print, 'There are',n_elements(goodnow),'OIII upper limit'
endif

;Overplot with Steidel et al, 2014 and Kewley et al, 2013 line. 
x=findgen(100)/50.-2.
ykewley = 0.61/(x+0.08)+1.1
ysteidel = 0.66/(x-0.31)+1.14
oplot,x,ykewley,color=60
oplot,x,ysteidel,color=90
stop

end
