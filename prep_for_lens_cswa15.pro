pro prep_for_lens_cswa15,type,name,path, Ha_file,oiii_file,metal_File
;chop input maps into many single frame maps whose coordinates match the reference image(fileref) e.g. hubble image used in Adi's lensing codes. Output is both tif and fits files with the same name with _interp at the end. 
  
; adapted from prep_for_lens.pro but this cswa15 doesn't have a ref file with RA-Dec coordinate.

  if type eq 'O3N2' then begin
     chop_map_15,path+Ha_file,1,name,'Ha' 
     chop_map_15,path+Ha_file,4,name,'Ha_detection' 
     chop_map_15,path+Ha_file,0,name,'kinematic' 
     chop_map_15,path+Ha_file,8,name,'kinematic_err' 
     chop_map_15,path+Ha_file,2,name,'veldisp' 
     chop_map_15,path+Ha_file,10,name,'veldisp_err' 
     chop_map_15,path+Ha_file,6,name,'NII' 
     chop_map_15,path+Ha_file,9,name,'Ha_err'
     chop_map_15,path+Ha_file,11,name,'NII_err'
     chop_map_15,path+Ha_file,12,name,'Ha_psf' 
     chop_map_15,path+OIII_file,1,name,'OIII' 
     chop_map_15,path+OIII_file,0,name,'kinematic_OIII' 
     chop_map_15,path+OIII_file,2,name,'veldisp_OIII' 
     chop_map_15,path+OIII_file,6,name,'Hb' 
     chop_map_15,path+OIII_file,9,name,'OIII_err'
     chop_map_15,path+OIII_file,11,name,'Hb_err'
     chop_map_15,path+metal_file,0,name,'N2index' 
     chop_map_15,path+metal_file,3,name,'N2index_err' 
     chop_map_15,path+metal_file,1,name,'O3N2index' 
     chop_map_15,path+metal_file,4,name,'O3N2index_err' 
     chop_map_15,path+metal_file,2,name,'Bayesianmetal' 
     chop_map_15,path+metal_file,5,name,'Bayesianmetal_err' 
     chop_map_15,path+metal_file,6,name,'metaltype' 
     
  endif

  if type eq 'N2' then begin
     chop_map_15,path+Ha_file,1,name,'Ha' 
     chop_map_15,path+Ha_file,4,name,'Ha_detection' 
     chop_map_15,path+Ha_file,0,name,'kinematic' 
     chop_map_15,path+Ha_file,8,name,'kinematic_err' 
     chop_map_15,path+Ha_file,2,name,'veldisp' 
     chop_map_15,path+Ha_file,10,name,'veldisp_err' 
     chop_map_15,path+Ha_file,6,name,'NII' 
     chop_map_15,path+Ha_file,12,name,'Ha_psf' 
     chop_map_15,path+metal_file,0,name,'N2index' 
     chop_map_15,path+metal_file,3,name,'N2index_err' 
     chop_map_15,path+metal_file,2,name,'Bayesianmetal' 
     chop_map_15,path+metal_file,5,name,'Bayesianmetal_err' 
     chop_map_15,path+metal_file,6,name,'N2metal_tag' 
  endif

 
end
