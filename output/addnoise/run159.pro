pro run159
for ii=0,10 do begin
filename='cswa159_Ha_tlc_Kc3_pipelinemosaic_sky_2hr_acube_noise'+strtrim(string(fix(ii)),2)+'.fits'
name = 'cswa159_'+strtrim(string(fix(ii)),2)
print, 'now doing file',name
;metalnew, type="N2",name=name,path='/scr2/nichal/workspace/output/addnoise/', Ha_file=filename,oiii_file='cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr_acube.fits',outpath='/scr2/nichal/workspace/output/addnoise/',/silent

metalfile = name+'_metallicitynew.fits'
;prep_for_lens,"N2",name,'/scr2/nichal/workspace/output/addnoise/', filename,'cswa159_Hb_tlc_Hn3_handmosaic_scaledsky_1hr_acube_aligned.fits',metalfile,'/scr2/nichal/workspace/imaging_data/cswa159/sdss_r.fits'

;;----------------------------
;;Below is from fixcswa159vel.pro
;;somehow the interpolation of the velocity map make an excessive rim
;files= file_Search('/scr2/nichal/workspace/output/forlensmodel/'+name+'/*_interp.fits')
;imref= readfits('/scr2/nichal/workspace/output/forlensmodel/CSWA159/CSWA159_kinematic_interp.fits')
;ind = array_indices(imref,findgen(n_elements(imref)))
;y = ind(1,*)
;bad1 = where(imref gt -35. and y ge 168)
;bad2 = where(abs(imref) le 12. and y lt 168)
;bad=[bad1,bad2]
;for n=0,n_Elements(files)-1 do begin
;   im = readfits(files(n),hdr)
;   im(bad) = 0.
;   writefits,files(n),im,hdr
;endfor
;;---------------------------------
;spawn,'\cp '+'/scr2/nichal/workspace/output/forlensmodel/'+name+'/*_interp.fits /scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/addnoise/'
;spawn,'\cp '+'/scr2/nichal/workspace/output/forlensmodel/'+name+'/*_interp.fits /scr2/nichal/workspace/massmodel/makemodel/BCK_scripts/LTMcswa159/outputmaps/addnoise/'

endfor

end
