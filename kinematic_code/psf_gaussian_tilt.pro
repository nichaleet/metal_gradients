function psf_gaussian_tilt,nx=nx,ny=ny,fwhm=fwhm,tilt=tilt,center=center
;fwhm = [a,b] a is the longer axis
;tilt in degrees
  X = FINDGEN(nx) # REPLICATE(1.0, ny)
  Y = REPLICATE(1.0, nx) # FINDGEN(ny)
; Define input function parameters:
  sigmax= fwhm(0)/2.355
  sigmay= fwhm(1)/2.355
  aAxis = sqrt(2.)*sigmax
  bAxis = sqrt(2.)*sigmay
  if keyword_set(center) then begin
     h=center(0)
     k=center(1)
  endif else begin
     h=round(nx/2.)
     k=round(ny/2.)
  endelse
  tilt = tilt*!PI/180  

; Create an ellipse:
  xprime = (X - h)*cos(tilt) - (Y - k)*sin(tilt)
  yprime = (X - h)*sin(tilt) + (Y - k)*cos(tilt)
  U = (xprime/aAxis)^2 + (yprime/bAxis)^2

; Create gaussian
  psf = exp(-U/2)/(2.*!pi*sigmax*sigmay)

; Chenck normallization
  tot = total(psf)
  psf = psf/tot

  return,psf
end
