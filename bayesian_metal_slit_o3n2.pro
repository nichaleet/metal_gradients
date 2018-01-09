Pro Bayesian_metal_slit_o3n2,name,n2_radius,n2_radius_err,o3hb_radius,o3hb_radius_err,radius,best_metal_arr,median_metal_arr,lowerbound_metal_arr,upperbound_metal_arr,metal_err_arr

;Calculated predicted OIII/Hb and NII/Ha from metallicity - Maiolino,2008
;INPUT: 
;if OIII image does not exist, input a 1./0. array with same size as Halpha array

;======================================================
set_plot,'x'

metal             = findgen(2300)/1000.+7. ; 12+log(O/H) = 7 to 9.3
x                 = metal-8.69
predictedN2  = -0.7732+(1.2357*x)-(0.2811*x^2)-(0.7201*x^3)-(0.3330*x^4)
predictedO3hb = 0.15490-(1.5031*x)-(0.9790*x^2)-(0.0297*x^3)  
;stop

;NII/Halpha
N2 = n2_radius
N2_err = n2_radius_err

;OIII/Hbeta
O3hb = o3hb_radius
O3hb_err = o3hb_radius_err

;Output arrays
best_metal_arr = fltarr(n_elements(radius))+1./0.  ;Least chisquare
median_metal_arr = best_metal_arr
lowerbound_metal_arr = best_metal_arr 
upperbound_metal_arr = best_metal_arr 

for ii=0,n_elements(radius)-1 do begin
   N2_current = N2(ii)
   N2_err_current = N2_err(ii)
   O3hb_current = O3hb(ii)
   O3hb_err_current = O3hb_err(ii) 
   ;Now  chisquare calculation.
   chisq = (N2_current-predictedN2)^2/N2_err_current^2+(O3hb_current-predictedO3hb)^2/O3hb_err_current^2 

   if where(finite(chisq) eq 0) ne [-1.] then stop
   print, 'min chisq is', min(chisq)
   print, '[NII]/Ha =',N2_current
   print, '[OIII]/Hb =',O3hb_current

   ConfidenceInterval,name,metal,chisq,'12+log(O/H)',216,x_peak,median_val,lower_val,upper_val
   best_metal_arr(ii) = x_peak(0)
   median_metal_arr(ii) = median_val 
   lowerbound_metal_arr(ii)= lower_val
   upperbound_metal_arr(ii) = upper_val
endfor
metal_err_arr = 0.5*(upperbound_metal_arr-lowerbound_metal_arr)
set_plot,'ps'
!p.multi = [0,1,1]
end
