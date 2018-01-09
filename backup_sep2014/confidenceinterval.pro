pro ConfidenceInterval,name,x_arr,chisq,x_paraname,type_current,x_peak,median_val,lower_val,upper_val
;calculate Posteria PDF and most likely values from chisquare
confidence = 0.68 ;1 standard deviation
  device, decomposed=0
;  loadct, 14
  !p.charsize = 1.0
  !p.color = 0 
  !p.background = 255
  !p.thick =2
  !x.thick =2
  !y.thick=2
  !p.multi=[0,2,1]

  
  ;change chisq to posteria probability distribution
  chisq = chisq-min(chisq)
  prob_x_arr    = exp(-0.5*chisq) ;probability (unscaled)
  area          = int_tabulated(x_arr,prob_x_arr,/double)
  prob_x_arr    = prob_x_arr/area    ;posteria pdf
  plot,x_arr,prob_x_arr,xtitle=x_paraname,ytitle='Posteria PDF',title=string(type_current)

  peaklocation = where(prob_x_arr eq max(prob_x_arr))
  x_peak = x_arr(peaklocation)
  print, 'most likely value is', x_peak
  cdf_arr = fltarr(n_elements(x_arr))
  for i =1,n_elements(x_arr)-1 do begin
     cdf_arr(i) = int_tabulated(x_arr[0:i],prob_x_arr[0:i])
  endfor
  plot,x_arr,cdf_arr,ytitle='cdf(x)',xtitle=x_paraname
  cdfpeak = cdf_arr(peaklocation)
;interpolation:
;1) median value
  median_pos = [max(where(cdf_arr-0.5 lt 0)),max(where(cdf_arr-0.5 lt 0))+1] ;array of two positions in x_arr where cdf is just below 0.5 and just above 0.5
  median_val = x_arr(median_pos(0))+(0.5-cdf_arr(median_pos(0)))*(x_arr(median_pos(1))-x_arr(median_pos(0)))/(cdf_arr(median_pos(1))-cdf_arr(median_pos(0)))

;2) upper bound
  ;max_prob = 0.5+confidence/2.
  ;if cdfpeak(0) gt 0.66 or cdfpeak(0) lt 0.34 then cdfpeak = 0.5
  max_prob = cdfpeak(0)+confidence/2.
  if max_prob lt 1. then begin
     upper_pos = [max(where(cdf_arr-max_prob lt 0)),max(where(cdf_arr-max_prob lt 0))+1] ;array of two positions in x_arr where cdf is just below confidence and just above confidence
     upper_val = x_arr(upper_pos(0))+(max_prob-cdf_arr(upper_pos(0)))*(x_arr(upper_pos(1))-x_arr(upper_pos(0)))/(cdf_arr(upper_pos(1))-cdf_arr(upper_pos(0)))
  endif else upper_val = x_arr(n_elements(x_arr)-1)
;3) lower bound
  ;min_prob = 0.5-confidence/2.
  min_prob = cdfpeak(0)-confidence/2.
  if min_prob gt 0. then begin
     lower_pos = [max(where(cdf_arr-min_prob lt 0)),max(where(cdf_arr-min_prob lt 0))+1] ;array of two positions in x_arr where cdf is just below confidence and just above confidence
     lower_val = x_arr(lower_pos(0))+(min_prob-cdf_arr(lower_pos(0)))*(x_arr(lower_pos(1))-x_arr(lower_pos(0)))/(cdf_arr(lower_pos(1))-cdf_arr(lower_pos(0)))
  endif else lower_val = x_arr(0)
  print, 'the confidence range of x is [',lower_val,',',upper_val,']'
  print, 'median of x value =', median_val
  ;wait, 0.5
end

