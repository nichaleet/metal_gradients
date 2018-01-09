pro unwrap,cube,filename

; Assume cube has dimensions lambda,x,y
cube = transpose(cube,[1,2,0])

sz = size(cube)
twod = fltarr(sz[3],sz[1]*sz[2])
k = 0d
for i = 0,sz[1]-1 do for j = 0,sz[2]-1 do begin
    twod[*,k] = cube[i,j,*]
    k = k+1
endfor

writefits,filename,twod

end
