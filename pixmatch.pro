pro pixmatch,pa_spec_1,scalex_1,scaley_1,raval_1,raVal_2,decval_1,decval_2,xref_1,yref_1,xref_2,yref_2,x_of_ref2_in1,y_of_ref2_in1


;This is for osiris data
;RA and DEC are in degrees
;Scalex_1 and scaley_1 are scales of image 1 in arcsec/pixel!!!

;;;;;;;;;;;;;;;;;;;
;Find pixel position of ra_2, dec_2 in image 1

GCIRC,1,RAval_1/15.,Decval_1,RAval_2/15.,Decval_2,ref_distance ; distance between the two reference points is in arcsec
POSANG,1,RAval_1/15.,Decval_1,RAval_2/15.,Decval_2,theta_prime ; Angle of the great circle containing ra_2,dec_2 from
;               the meridian containing [ra_1, dc_1], in the sense north
;               through east rotating about [ra1, dc1].

theta = theta_prime-PA_spec_1  ; in degrees
x_of_ref2_in1  = ref_distance*cos(theta*!dpi/180.)/scalex_1+xref_1
y_of_ref2_in1  = ref_distance*sin(theta*!dpi/180.)/scaley_1+yref_1

print,'(x,y) in image2:', xref_2, yref_2
print,'matched(x,y) in image1', x_of_ref2_in1,y_of_ref2_in1

x_of_ref2_in1  = round(x_of_ref2_in1)
y_of_ref2_in1  = round(y_of_ref2_in1)

;stop
end
