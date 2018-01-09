Pro Bayesian_metal,name,OIII,hbeta,NII,halpha,OIII_err,hbeta_err,NII_err,halpha_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian,silent=silent

;Calculated predicted OIII/Hb and NII/Ha from metallicity - Maiolino,2008
;INPUT: 
;if OIII image does not exist, input a 1./0. array with same size as Halpha array

;======================================================
set_plot,'x'
metal             = findgen(2300)/1000.+7. ; 12+log(O/H) = 7 to 9.3
x                 = metal-8.69
predictedlog_NII_Ha  = -0.7732+(1.2357*x)-(0.2811*x^2)-(0.7201*x^3)-(0.3330*x^4)
predictedlog_OIII_Hb = 0.15490-(1.5031*x)-(0.9790*x^2)-(0.0297*x^3)  
predictedNII_Ha = 10.^predictedlog_NII_Ha
predictedOIII_Hb = 10.^predictedlog_OIII_Hb
;stop
;nan where oiii/hbeta is negative. This will give more importance to NII/Halpha since if it's negative then the only result given is the highest metallicity. The predicted OIII/Hb when it turns over to equal to the original point at minimum z is ~1.1. The min(oiii_hb(where(oiii_hb gt 0.))) ~ 1.4. Hence there is no way that the non-Hb detection is gonna be because of high metallicity.
bad=where(hbeta lt 0.)
oiii(bad) = 1./0.
hbeta(bad) = 1./0.
oiii_err(bad) = 1./0.
hbeta_err(bad) = 1./0.
if ~keyword_set(silent) then print,'total of ',n_elements(bad),' oiii/hb pixels had been changed'
;stop
;find detected NII pixels
NII_current = NII
NII_upper = NII_err+NII
change2upper = where(NII lt 0. and NII_upper gt 0.)
sizeHa = size(halpha)
upper_tags = fltarr(sizeHa(1),sizeHa(2))
upper_tags(where(finite(NII_current))) = 8.
if change2upper(0) ne -1 then begin
   ;NII_current(change2upper) = NII_upper(change2upper)
   upper_tags(change2upper) = 4.
endif
;stop
nan = where(NII_upper lt 0.)
upper_tags(nan) = 2.
;NII_current(nan) = 1./0.

;Now tag=8 is detection, tag =4 is upper limit detection, 
                                ;tag=2 is non detection
detectedNIIpix = where(upper_tags eq 8. or upper_tags eq 4.)

;find detected Hbeta pixels
Hbeta_current = Hbeta
Hbeta_upper = Hbeta_err+Hbeta
change2upper = where(Hbeta lt 0. and Hbeta_upper gt 0.)
HbetaUpper_tags = fltarr(sizeHa(1),sizeHa(2)) 
finitepix = where(finite(hbeta))     
if finitepix(0) ne -1 then HbetaUpper_tags(finitepix) = 27.

if change2upper(0) ne -1 then begin
   ;Hbeta_current(change2upper) = Hbeta_upper(change2upper)
   HbetaUpper_tags(change2upper) = 9.
endif
nan = where(Hbeta_upper lt 0.)
if nan(0) ne -1 then HbetaUpper_tags(nan) = 3.
;Hbeta_current(nan) = 1./0.
;Now tag=27 is detection, tag =9 is upper limit detection, 
                                ;tag=3 is non detection
detectedHbetapix = where(HbetaUpper_tags eq 27. or HbetaUpper_tags eq 9.)
if detectedHbetapix(0) eq -1 then detectedHbetapix = []

;NII/Halpha
NII_Ha = NII/halpha
NII_Ha_err = sqrt((NII_err/NII)^2+(Halpha_err/Halpha)^2)
;sigma(logA) = sigma(a)/A/log10
;OIII/Hbeta
OIII_Hb = OIII/Hbeta
OIII_Hb_err = sqrt((OIII_err/OIII)^2+(Hbeta_err/Hbeta)^2)

;Output arrays
best_metal_arr = halpha*1/0.  ;Least chisquare
median_metal_arr = halpha*1/0.   ;CDF = 0.5
lowerbound_metal_arr = halpha*1/0.
upperbound_metal_arr = halpha*1/0.
typebayesian = halpha*1./0.

;Find pixels to run for loop for. (Pixels with either Hb or NII detections)
goodpix = SetUnion(detectedNIIpix,detectedHbetapix)
loadct, 14
;for loop over all NII or Hb detection pixels
for ii=0,n_elements(goodpix)-1 do begin
   NII_Ha_current = NII_Ha(goodpix(ii))
   NII_Ha_err_current = NII_Ha_err(goodpix(ii))
   OIII_Hb_current = OIII_Hb(goodpix(ii))
   OIII_Hb_err_current = OIII_Hb_err(goodpix(ii)) 
   NII_current = NII(goodpix(ii))
   Hb_current  = HBeta(goodpix(ii))
                                ;sometimes the errors are -Nan instead of Inf. So, we have to change them to infinity
   if finite(NII_Ha_err_current) eq 0. then NII_Ha_err_current = 1./0.
   if finite(OIII_Hb_err_current) eq 0. then OIII_Hb_err_current = 1./0.

                                ;Now separating chisquare calculation. Both Hb and NII, NII only, Hb only
   typeofchi = [finite(NII_current),finite(Hb_current)] 
   if array_equal(typeofchi,[1,1]) then begin
      chisq = (NII_Ha_current-predictedNII_Ha)^2/NII_Ha_err_current^2+(OIII_Hb_current-predictedOIII_Hb)^2/OIII_Hb_err_current^2 
   endif
   if array_equal(typeofchi,[0,1]) then chisq = (OIII_Hb_current-predictedOIII_Hb)^2/OIII_Hb_err_current^2
      
   if array_equal(typeofchi,[1,0]) then chisq = (NII_Ha_current-predictedNII_Ha)^2/NII_Ha_err_current^2
   if array_equal(typeofchi,[0,0]) then begin 
      if  ~keyword_set(silent) then print,'Both Ha and OIII are not detected in this pixel. Something is wrong.'
      stop
   endif

   if where(finite(chisq) eq 0) ne [-1.] then stop
   if  ~keyword_set(silent) then begin
   print, 'min chisq is', min(chisq)
   print, '[NII]/Ha =',NII_Ha_current
   print, '[OIII]/Hb =',OIII_Hb_current
   endif
   type_current = Upper_tags(goodpix(ii))*Hbetaupper_tags(goodpix(ii))
   if type_current eq 0. then type_current = Upper_tags(goodpix(ii))+Hbetaupper_tags(goodpix(ii))

   typeBayesian(goodpix(ii)) = type_current
   if  ~keyword_set(silent) then print, 'type of chisq is', type_current
   ConfidenceInterval,name,metal,chisq,'12+log(O/H)',type_current,x_peak,median_val,lower_val,upper_val,silent=silent
   ;if type_current ne 8. then stop
   best_metal_arr(goodpix(ii)) = x_peak(0)
   median_metal_arr(goodpix(ii)) = median_val 
   lowerbound_metal_arr(goodpix(ii))= lower_val
   upperbound_metal_arr(goodpix(ii)) = upper_val
   if type_current eq 72 then stop
endfor
;stop
set_plot,'ps'
!p.multi = [0,1,1]
end
