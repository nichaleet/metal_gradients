function makeborder,cube,bordersize
  addpixel = 2.*bordersize
  sizecube = size(cube,/dimensions)
  newcube = fltarr(sizecube[0]+addpixel,sizecube[1]+addpixel,sizecube[2])+1./0.
  newcube[bordersize:sizecube[0]+bordersize-1,bordersize:sizecube[1]+bordersize-1,*]=cube
  return,newcube
end


pro splitfile

;File 1
files=file_search('*ratios.fits',/FOLD_CASE)
for i=0,n_elements(files)-1 do begin
   name = strmid(files(i),0,strpos(files(i),'_'))
   print,name
   cube=readfits(files(i))
   
   if name eq '1148' then begin
      cube[64,75,*]=1./0.
      cube[65,9,*]=1./0.
   endif
   if name eq '0744' then begin
      cube = makeborder(cube,20.)
      sizecube = size(cube,/dimensions)
      badpix = where(cube[*,*,0] lt -200.)
      badpixarr = fltarr(sizecube(0),sizecube(1))+1.
      badpixarr(badpix) = 0.
      badpixnow = rebin(badpixarr,sizecube(0),sizecube(1),sizecube(2))
      cube = cube/badpixnow
   endif

   if name eq '1038' then cube = makeborder(cube,20.)
   writefits,name+'_vel.fits',cube[*,*,0]
   writefits,name+'_Ha.fits',cube[*,*,1]
   writefits,name+'_NII_Ha.fits',cube[*,*,2]
   writefits,name+'_OIII_Hb.fits',cube[*,*,4]

;File 2
   name_errfile = 'acube_'+name+'_ha_error.fits'
   cube2 = readfits(name_errfile)
   if name eq '0744' then cube2 = makeborder(cube2,20.)
   if name eq '1038' then cube2 = makeborder(cube2,20.)

   writefits,name+'_velerr.fits',cube2[*,*,0]
   writefits,name+'_Haerr.fits',cube2[*,*,1]
   writefits,name+'_veldisperr.fits',cube2[*,*,2]

; calculate NII, NII err, NII_Ha err
   NII_Ha = cube[*,*,2]
   Ha = cube[*,*,1]
   Ha_err = cube2[*,*,1]
   NII = NII_Ha*Ha
   NII_err = Ha_err/Ha*NII  ;Ha_err
   NII_Ha_err = NII_Ha*sqrt(2.*(Ha_err/Ha)^2)
   writefits,name+'_NII.fits',NII
   writefits,name+'_NII_err.fits',NII_err
   writefits,name+'_NII_Ha_err.fits',NII_Ha_err

;File 3
   name_hafile = 'acube_'+name+'_ha.fits'
   cube3 = readfits(name_hafile)
   if name eq '0744' then cube3 = makeborder(cube3,20.)
   if name eq '1038' then cube3 = makeborder(cube3,20.)

   writefits,name+'_hadetect.fits',cube3[*,*,4]
   writefits,name+'_veldisp.fits',cube3[*,*,2]
endfor
stop
end 
