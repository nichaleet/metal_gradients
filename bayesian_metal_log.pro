Pro Bayesian_metal_log,name,o3hb,n2,o3hb_err,n2_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian,silent=silent

;Calculated predicted OIII/Hb and NII/Ha from metallicity - Maiolino,2008
;INPUT: 
;if OIII image does not exist, input a 1./0. array with same size as Halpha array

;======================================================

;check for negative values
if (where(n2 lt 0.))[0] ne -1 then begin
	n2(where(n2 lt 0.)) = 1./0.
	n2_err(Where(n2 lt 0.)) = 1./0.
endif
if (where(o3hb lt 0.))[0] ne -1 then begin
	o3hb(where(o3hb lt 0.)) = 1./0.
	o3hb_err(where(o3hb lt 0.)) = 1./0.
endif 

set_plot,'x'
metal             = findgen(2300)/1000.+7. ; 12+log(O/H) = 7 to 9.3
x                 = metal-8.69
predictedlog_NII_Ha  = -0.7732+(1.2357*x)-(0.2811*x^2)-(0.7201*x^3)-(0.3330*x^4)
predictedlog_OIII_Hb = 0.15490-(1.5031*x)-(0.9790*x^2)-(0.0297*x^3)  

detectedNIIpix = where(finite(n2)) 
detectedHbetapix = where(finite(o3hb)
if detectedHbetapix(0) eq -1 then detectedHbetapix = []

;Output arrays
best_metal_arr = n2*1/0.  ;Least chisquare
median_metal_arr = n2*1/0.   ;CDF = 0.5
lowerbound_metal_arr = n2*1/0.
upperbound_metal_arr = n2*1/0.
typebayesian = n2*1./0.

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
      if not (silent) then print,'Both Ha and OIII are not detected in this pixel. Something is wrong.'
      stop
   endif

   if where(finite(chisq) eq 0) ne [-1.] then stop
   if not (silent) then begin
   print, 'min chisq is', min(chisq)
   print, '[NII]/Ha =',NII_Ha_current
   print, '[OIII]/Hb =',OIII_Hb_current
   endif
   type_current = Upper_tags(goodpix(ii))*Hbetaupper_tags(goodpix(ii))
   if type_current eq 0. then type_current = Upper_tags(goodpix(ii))+Hbetaupper_tags(goodpix(ii))

   typeBayesian(goodpix(ii)) = type_current
   if not (silent) then print, 'type of chisq is', type_current
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
