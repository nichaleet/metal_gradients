pro prep_for_lens,type,name,path, Ha_file,oiii_file,metal_File,fileref
;chop input maps into many single frame maps whose coordinates match the reference image(fileref) e.g. hubble image used in Adi's lensing codes. Output is both tif and fits files with the same name with _interp at the end. 


  if type eq 'O3N2' then begin
     chop_map,path+Ha_file,1,name,'Ha',fileref
     chop_map,path+Ha_file,0,name,'kinematic',fileref
     chop_map,path+Ha_file,2,name,'veldisp',fileref
     chop_map,path+Ha_file,6,name,'NII',fileref
     chop_map,path+OIII_file,1,name,'OIII',fileref
     chop_map,path+OIII_file,0,name,'kinematic_OIII',fileref
     chop_map,path+OIII_file,2,name,'veldisp_OIII',fileref
     chop_map,path+OIII_file,6,name,'Hb',fileref
     chop_map,path+metal_file,0,name,'N2metal',fileref
     chop_map,path+metal_file,1,name,'O3N2metal',fileref
     chop_map,path+metal_file,2,name,'Bayesianmetal',fileref
     chop_map,path+metal_file,6,name,'metaltype',fileref
     
  endif

  if type eq 'N2' then begin
     chop_map,path+Ha_file,1,name,'Ha',fileref
     chop_map,path+Ha_file,0,name,'kinematic',fileref
     chop_map,path+Ha_file,2,name,'veldisp',fileref
     chop_map,path+Ha_file,6,name,'NII',fileref
     chop_map,path+metal_file,0,name,'N2metal',fileref
     chop_map,path+metal_file,1,name,'Bayesianmetal',fileref
     chop_map,path+metal_file,4,name,'N2metal_tag',fileref
  endif
end
