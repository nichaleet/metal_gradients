pro tiftofits,fileref

sourcetiffiles = file_search('source*.tif')
outname = repstr(sourcetiffiles,'tif','fits')
print, sourcetiffiles
ref = readfits(fileref,header)

restore,'fitstifconver.sav'  ;This will give tiffiles and linvararr (intercept and slope)
print, 'Make sure the list belows are the same!!!!'
print, sourcetiffiles
print, tiffiles
stop
;sourceHa.tif sourceN2metal.tif sourceO3N2metal.tif sourceOrig.tif
;sourcekinematic.tif sourceveldisp.tif
slope = linvararr(1,*);[0.000622361,0.0348441,0.0327519,1.,1.80299,1.76104]
intercept = linvararr(0,*);[-2.31040e-8,-3.24031e-6,-6.97351e-7,0.,-169.481,-4.97651e-5]

for ii=0, n_elements(sourcetiffiles)-1 do begin

img = read_tiff(sourcetiffiles(ii))
;if ii eq 4 then background=where(img eq 94 or img eq 0)
img = img*slope(ii)+intercept(ii)
imgcopy = img
mode = mode(imgcopy)
background = where(img eq mode)
print, mode
img(background) = 1./0.
plot,img,psym=1
wait,2
;if tiffiles(ii) eq 'sourceOrig.tif' then img = transpose(img,[1,2,0])

writefits,outname(ii),img,header
endfor

end
