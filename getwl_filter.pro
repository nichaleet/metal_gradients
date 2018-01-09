function getwl_filter,filter

if filter eq 'Jbb' then wlarr =1.180+1.5D-4*dindgen(1574)else $
if filter eq 'Jn1' then wlarr =1.174+1.4987D-4*dindgen(388)else $
if filter eq 'Jn2' then wlarr =1.228+1.5D-4*dindgen(408)else $
if filter eq 'Hbb' then wlarr =1.473+2.D-4*dindgen(1651)else $
if filter eq 'Hn1' then wlarr =1.466+2.D-4*dindgen(376)else $
if filter eq 'Hn2' then wlarr =1.532+2.D-4*dindgen(391)else $
if filter eq 'Hn3' then wlarr =1.594+2.D-4*dindgen(411)else $
if filter eq 'Hn4' then wlarr =1.652+2.D-4*dindgen(426)else $
if filter eq 'Kn1' then wlarr =1.955+2.5D-4*dindgen(401)else $
if filter eq 'Kn2' then wlarr =2.036+2.5D-4*dindgen(421)else $
if filter eq 'Kc5' then wlarr =2.292+2.5D-4*dindgen(465)else $
if filter eq 'Kc3' then wlarr =2.121+2.5D-4*dindgen(433)else begin 
   print, 'no filter match'   
   stop
endelse

return,wlarr

end
