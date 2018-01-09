Pro Bayesian_metal_log,name,o3hb,o3hb_err,n2,n2_err,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,typebayesian,silent=silent

;Calculated predicted OIII/Hb and NII/Ha from metallicity - Maiolino,2008
;INPUT: 
;if OIII image does not exist, input a 1./0. array with same size as Halpha array

;======================================================
set_plot,'x'
metal             = findgen(2300)/1000.+7. ; 12+log(O/H) = 7 to 9.3
x                 = metal-8.69
predictedn2  = -0.7732+(1.2357*x)-(0.2811*x^2)-(0.7201*x^3)-(0.3330*x^4)
predictedO3Hb = 0.15490-(1.5031*x)-(0.9790*x^2)-(0.0297*x^3)  
if n_elements(silent) eq 0 then silent = 0

detectedNIIpix = where(finite(n2))
if detectedniipix(0) eq -1 then detectedniipix=[]
detectedo3hbpix = where(finite(o3hb))
if detectedo3hbpix(0) eq -1 then detectedo3hbpix = []

;Output arrays
best_metal_arr = n2*1/0.  ;Least chisquare
median_metal_arr = n2*1/0.   ;CDF = 0.5
lowerbound_metal_arr = n2*1/0.
upperbound_metal_arr = n2*1/0.
typebayesian = n2*1./0.

;Find pixels to run for loop for. (Pixels with either Hb or NII detections)
goodpix = SetUnion(detectedNIIpix,detectedo3hbpix)
loadct, 14
;for loop over all NII or Hb detection pixels
for ii=0,n_elements(goodpix)-1 do begin
   N2_current = N2(goodpix(ii))
   N2_err_current = N2_err(goodpix(ii))
   O3hb_current = O3hb(goodpix(ii))
   O3hb_err_current = O3hb_err(goodpix(ii)) 
   ;sometimes the errors are -Nan instead of Inf. So, we have to change thm to infinity
   if finite(N2_err_current) eq 0. then N2_err_current = 1./0.
   if finite(O3hb_err_current) eq 0. then O3hb_err_current = 1./0.

   ;Now separating chisquare calculation. Both Hb and NII, NII only, Hb only
   typeofchi = [finite(N2_current),finite(O3Hb_current)] 
   if array_equal(typeofchi,[1,1]) then begin
      chisq = (N2_current-predictedN2)^2/N2_err_current^2+(O3hb_current-predictedO3hb)^2/O3hb_err_current^2 
   endif
   if array_equal(typeofchi,[0,1]) then chisq = (O3hb_current-predictedO3hb)^2/O3hb_err_current^2
      
   if array_equal(typeofchi,[1,0]) then chisq = (N2_current-predictedN2)^2/N2_err_current^2
   if array_equal(typeofchi,[0,0]) then begin 
      if not (silent) then print,'Both Ha and OIII are not detected in this pixel. Something is wrong.'
      stop
   endif

   if where(finite(chisq) eq 0) ne [-1.] then stop
   if not (silent) then begin
      print, 'min chisq is', min(chisq)
      print, '[NII]/Ha =',N2_current
      print, '[OIII]/Hb =',O3hb_current
   endif

   typecurrent = total([2.*finite(n2_current),3.*finite(o3hb_current)])
   ConfidenceInterval,name,metal,chisq,'12+log(O/H)',typecurrent,x_peak,median_val,lower_val,upper_val
   best_metal_arr(goodpix(ii)) = x_peak(0)
   median_metal_arr(goodpix(ii)) = median_val 
   lowerbound_metal_arr(goodpix(ii))= lower_val
   upperbound_metal_arr(goodpix(ii)) = upper_val
endfor
;stop
set_plot,'ps'
!p.multi = [0,1,1]
end
