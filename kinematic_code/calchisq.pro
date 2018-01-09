function calchisq,ymodel,y,yerr,nparam
;return reduced chisq
if n_Elements(ymodel) ne n_elements(y) then begin
   print, 'Number of elements are not equal'
   stop
endif

chisq = total((ymodel-y)^2/yerr^2)
chisq = chisq/(n_elements(ymodel)-nparam-1.)
;stop
return,chisq

end
