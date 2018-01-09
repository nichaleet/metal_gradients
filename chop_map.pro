pro chop_map,file,framenumber,name,framename,fileref

incube = readfits(file, header)
sizecube = size(incube)
if framenumber lt sizecube(3) then begin
   image = incube(*,*,framenumber)   
   image(where(finite(image) eq 0.)) = 0.
   if name eq 'CSWA11_firstgal' and framename eq 'N2index' then image(where(abs(image) gt 1.)) = 0.
   if name eq 'CSWA11_secondgal' and framename eq 'N2index' then image(where(abs(image) gt 1.)) = 0.
   if name eq 'CSWA11_firstgal' and framename eq 'N2index_err' then image(where(abs(image) gt 1.)) = 0.   
   if name eq 'CSWA11_secondgal' and framename eq 'N2index_err' then image(where(abs(image) gt 1.)) = 0.   
   if name eq 'CSWA159' and framename eq 'kinematic' then image[35:36,41]=0. 
   dir_name = '/scr2/nichal/workspace/output/forlensmodel/'+name
   file_mkdir,dir_name
   filenameout = dir_name+'/'+name+'_'+framename+'.fits'
   writefits,filenameout,image,header

   interp_map,filenameout,fileref
endif else print, 'frame', framenumber,' does not exist.'
end
