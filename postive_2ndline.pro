pro postive_2ndline
window, 0 , xsize =700,ysize=700.
!p.multi = [0,1,2]
restore,'positive_2ndline.sav' ;fitnegativestr
restore,'fitstruct.sav'   ;fitstruct

for i=0, n_elements(fitpositivestr)-1 do begin
NIInow = fitpositivestr(i)
xcoord = NIInow.x
ycoord = NIInow.y
sigma  = NIInow.sigma
if ycoord lt 20 then goto, skip

matchHA = where(fitstruct.x eq xcoord and fitstruct.y eq ycoord)
if matchHa(0) eq -1  then stop else begin
   if n_Elements(matchHA) gt 1. then matchHA = matchHA(n_elements(matchHA)-1)
   Hanow = fitstruct(matchHA)
   sigHA = Hanow.sigma
   plot, Hanow.wl,Hanow.data,psym=1,title=string(xcoord)+string(ycoord)+string(sigHa)
   oplot,Hanow.wl,Hanow.fit,color=255
   redshift = Hanow.wl(where(Hanow.fit eq max(Hanow.fit)))/.65628-1.

   plot, NIInow.wl,NIInow.data,psym=1,title=string(xcoord)+string(ycoord)+string(sigma)
   oplot,NIInow.wl,NIInow.fit,color=255
   ;oplot,[expectedNII,expectedNII],!y.crange,linestyle=2
   oplot,[1.6752,1.6752],!y.crange,linestyle=2
   oplot,[1.6758,1.6758],!y.crange,linestyle=2
   
endelse
stop
skip: continue
endfor


end
