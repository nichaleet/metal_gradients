pro make_hst_small_fit

;;a773
;image = readfits('a773.fits',header,exten_no=1)
;hextract,image,header,newim,newhd,832,872,622,686

;;copy the region where osiris frame is
;image_osiris = image(832:872,622:686)

;;set the blank area to nan
;newim(0:10,*) = 1./0.
;newim(*,0:35) = 1./0.
;newim(23:40,*) = 1./0.
;newim(*,55:64) = 1./0.

;writefits,'a773_small.fits',newim,newhd

;;macs1133
;image = readfits('MACSJ1133.2+5008_F606W.fits',header,exten_no=1)
;hextract,image,header,newim,newhd,1908,2100,2258,2440


;;INSERTING A SQUARE OF OSIRIS image in the Hubble small file
;array of x,y for the newim
;ind = array_indices(newim,where(newim))
;xind = ind(0,*)
;yind = ind(1,*)

;line1 = where(abs((xind-66.)*49./33.-yind+37.) le 2. and xind ge 66. and xind le 99.) 
;line2 = where(abs((xind-181.)*(-59./82.)-yind+27.) le 2. and xind ge 99. and xind le 181.)
;line3 = where(abs((xind-114.)*(-37./48.)-yind) le 2. and xind ge 66. and xind le 114.)  
;line4 = where(abs((xind-162.)*27./19.-yind) le 2. and xind ge 162. and xind le 181.) 

;square = [line1,line2,line3,line4]
;newim(square)=0.1
;writefits,'macs1133_small.fits',newim,newhd

;extract osiris frame (kind of)
;hextract,image,header,newim,newhd,1988,2086,2260,2330
;writefits,'macs1133_small_osiris.fits',newim,newhd


image = readfits('hlsp_clash_hst_acs_macs0717_f814w_v1_drz.fits',header);,exten_no=1)
hextract,image,header,newim,newhd,2483,2544,2511,2571

;;set the blank area to nan
;newim(0:10,*) = 1./0.
;newim(*,0:35) = 1./0.
;newim(23:40,*) = 1./0.
;newim(*,55:64) = 1./0.

writefits,'macs0717_small.fits',newim,newhd

stop
end
