pro a773_makemap
; xsky [ysky] is list of x [y] values of points on sky - these transform
; to position xsource [ysource] when the lensing potential is taken
; into account.
readcol, 'A773_source.dat',ids,xsource,ysource ;(id, RA,DEC)
readcol, 'A773_Osiris.mul',idi,xsky,ysky       ;image plane (idi, RA, DEC)

imref=readfits('/scr2/nichal/workspace/output/abell773_Ha_tlc_Kc3_handmosaic_sky_330hr_acube.fits',hdr)
imref=imref[*,*,1]
extast,hdr,astr
ad2xy,xsky,ysky,astr,x,y

sizeim = size(imref,/dimensions)
newmapx = fltarr(sizeim(0),sizeim(1),10)+0.
newmapy = newmapx
countmap = fltarr(sizeim)
for i=0,n_elements(x)-1 do begin
   currentx = fix(x(i))
   currenty = fix(y(i))
   if currentx ge sizeim(0) or currentx lt 0 then goto, skip
   if currenty ge sizeim(1) or currenty lt 0 then goto, skip
   num_in_pix = countmap[currentx,currenty]
   newmapx[currentx,currenty,num_in_pix] = xsource(i)
   newmapy[currentx,currenty,num_in_pix] = ysource(i)
   countmap[currentx,currenty]+=1   
   if num_in_pix gt 9 then stop
   skip:continue
endfor
stop
   newmapx=total(newmapx,3)/countmap
   newmapy=total(newmapy,3)/countmap

finalmap=[[[imref]],[[newmapx]],[[newmapy]]]
writefits,'a773_model.fits',finalmap,hdr
stop
end
