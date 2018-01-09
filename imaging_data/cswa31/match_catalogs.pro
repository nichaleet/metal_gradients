;This pro match the cswa31 catalogs based on the r band image. It'll check the position of objects in r catalog and check the check.fits image of other bands for the id

pro match_catalogs

gcat = 'cswa31_g_noheader.cat'
icat = 'cswa31_i_noheader.cat'
rcat = 'cswa31_r_noheader.cat'

checki = readfits('checki.fits')
checkg = readfits('checkg.fits')

readcol,rcat,r1,R2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28


readcol,gcat,g1,G2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g25,g26,g27,g28

gold = [[g1],[G2],[g3],[g4],[g5],[g6],[g7],[g8],[g9],[g10],[g11],[g12],[g13],[g14],[g15],[g16],[g17],[g18],[g19],[g20],[g21],[g22],[g23],[g24],[g25],[g26],[g27],[g28]]
gold = transpose(gold)

readcol,icat,i1,I2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28

iold = [[i1],[I2],[i3],[i4],[i5],[i6],[i7],[i8],[i9],[i10],[i11],[i12],[i13],[i14],[i15],[i16],[i17],[i18],[i19],[i20],[i21],[i22],[i23],[i24],[i25],[i26],[i27],[i28]]
iold = transpose(iold)

nobj = n_elements(r1)
inew = fltarr(28,nobj)
gnew = fltarr(28,nobj)

for ii=0,nobj-1 do begin
   id = r1(ii) 
   x  = round(r2(ii))
   y  = round(r3(ii))
;G image
   g_id = checkg(x,y)
   if g_id ne 0 then begin  ;0 is no match found
      g_pos = where(g1 eq g_id)
      gnew(*,ii) = gold(*,g_pos)
   endif else  for jj=0,27 do gnew(jj,ii) = 999
      

;I image
   i_id = checki(x,y)
   if i_id ne 0 then begin  ;0 is no match found
      i_pos = where(i1 eq i_id)
      inew(*,ii) = iold(*,i_pos)
   endif else for jj=0,27 do inew(jj,ii) = 999

                                ;if id eq 193 then stop
endfor

openw,1,'cswa31_g_noheader_match.cat'
printf,1,gnew,FORMAT='(I5,F11.4,F11.4,F12.7,F12.7,F9.4,F9.4,F9.4,F9.4,F9.4,F9.4,F18.3,F10.4,F18.3,F10.5,F18.3,F10.4,F8.2,F8.3,F8.2,F8.3,F8.3,F8.2,F8.3,F8.3,I5,F9.4,F10.4)'    
close,1

openw,1,'cswa31_i_noheader_match.cat'
printf,1,inew,FORMAT='(I5,F11.4,F11.4,F12.7,F12.7,F9.4,F9.4,F9.4,F9.4,F9.4,F9.4,F18.3,F10.4,F18.3,F10.5,F18.3,F10.4,F8.2,F8.3,F8.2,F8.3,F8.3,F8.2,F8.3,F8.3,I5,F9.4,F10.4)'   
close,1


end
