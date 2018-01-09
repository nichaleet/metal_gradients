pro tiftofits,fileref

sourcetiffiles= ['sourceBayesianmetal.tif','sourceha.tif','sourcehb.tif','sourcen2index.tif','sourcenii.tif','sourceo3n2index.tif','sourceoiii.tif','sourcekinematic.tif','sourcemetaltype.tif','sourceveldisp.tif']
outname = repstr(sourcetiffiles,'tif','fits')
print, sourcetiffiles
ref = readfits(fileref,header)


restore,'fitstifconver.sav'  ;This will give tiffiles and linvararr (intercept and slope)
print, 'Make sure the lists below are the same!!!!'
for ii=0,n_elements(sourcetiffiles)-1 do print, sourcetiffiles(ii),' ',tiffiles(ii)
stop
;sourceHa.tif sourceN2metal.tif sourceO3N2metal.tif sourceOrig.tif
;sourcekinematic.tif sourceveldisp.tif
slope = linvararr(1,*)
intercept = linvararr(0,*)


for ii=0, n_elements(sourcetiffiles)-1 do begin

img = read_tiff(sourcetiffiles(ii))

;if sourcetiffiles(ii) eq 'sourceNII.tif' then img(where(img eq 0.)) = 148.
;if sourcetiffiles(ii) eq 'sourcekinematic.tif' then img(where(img eq 0.)) = 135.
img = img*slope(ii)+intercept(ii)


imgcopy = img
mode = mode(imgcopy)
background = where(img eq mode)
print,sourcetiffiles(ii), mode,n_elements(background)
img(background) = 1./0.
print, 'background value', mode
plot,img,psym=1
wait,2
;if tiffiles(ii) eq 'sourceOrig.tif' then img = transpose(img,[1,2,0])

writefits,outname(ii),img,header
endfor

end
