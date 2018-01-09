pro chop_map_15,file,framenumber,name,framename

incube = readfits(file, header)
image = incube(*,*,framenumber)
image(where(finite(image) eq 0.)) = 0.

dir_name = '/scr2/nichal/workspace/output/forlensmodel/'+name
file_mkdir,dir_name
filenameout = dir_name+'/'+name+'_'+framename+'.fits'
writefits,filenameout,image,header
interp_map_15,filenameout

end
