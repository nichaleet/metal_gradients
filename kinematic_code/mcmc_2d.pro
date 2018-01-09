pro mcmc_2d,x1,x2,y,yerr,nparam,nstep,p0,stepsizes,rangemin,rangemax,paraname,kernel,ind,bad,name,pixscale,fixparam,returnvalues  
   ;The input is x1, x2 and y.
   ;Fit y with input data x1 and x2 with a function of n parameters.
   ;Return values is an array of size 3*nparam contain 
    ;lower limit,best estimate, and upper limit values for each parameter.
   ;Functions that are needed in this code are calymodel_2d_conv and calchisq

ngoodparam  = n_elements(where(fixparam ne 0.))
ymodel0 = calymodel_2d_conv(x1,x2,p0,kernel,ind,bad,0,name)
oldchisq = calchisq(ymodel0,y,yerr,ngoodparam)
maxradius = calmaxradius(x1,x2,p0,kernel,ind,bad,0,pixscale)
if maxradius gt 20. then stop,'ERROR You picked too low inclination'

print,'Chisq of the initial parameters:', oldchisq

randomparam = fix(randomu(seed,nstep)*ngoodparam) 
goodparam  = where(fixparam ne 0.)
     ;choose a parameter to perturbed in each step 
randomstepsize = randomn(seed,nstep) 
     ;randomn returns a random numbers from gaussian dist 
     ;with mean of zero and sigma=1.1
randomaccept = randomu(seed,nstep)

parr = fltarr(nparam,nstep+1) ;keep parameters for each step
chisqarr = fltarr(nstep+1)    ;keep chisquare for each step
pold = p0
parr(*,0) = p0
pcurrent = p0
chisqarr(0) = oldchisq
naccepts = 0UL  ;To count how many steps have been accepted
nrejects = 0UL  ;To count how many steps have been rejected
     ;(note) we want about 20%-50% acceptance rate

for ii=0UL,nstep-1 do begin
   if (ii mod 10000.) eq 0. then begin
      print, 'Doing step',ii
      print, 'i,theta,xc,yc,rt,vc,v0 =', pcurrent
      print,'chisq =' ,oldchisq
   endif
   chosenparam = goodparam(randomparam(ii)) ;randomly chosen parameter to be perturbed
   pcurrent(chosenparam) = pcurrent(chosenparam)+stepsizes(chosenparam)*randomstepsize(ii)
   ymodel = calymodel_2d_conv(x1,x2,pcurrent,kernel,ind,bad,ii,name)                      
      ;since P(m|D) = P(D|m)*P(m|I). 
      ;For the model that is out of range, we plus infinity to the chisq
   newchisq = calchisq(ymodel,y,yerr,ngoodparam)
   if pcurrent(chosenparam) gt rangemax(chosenparam) or pcurrent(chosenparam) lt rangemin(chosenparam) then newchisq=1.e8
                                ;another constraint if max(radius) of detected pixel is greater than 20 kpc then the model is bad
   maxradius = calmaxradius(x1,x2,pcurrent,kernel,ind,bad,ii,pixscale)
   if maxradius gt 20. then newchisq = 1.e8

;   print, pcurrent, newchisq   
   if newchisq lt oldchisq then begin 
      naccepts = naccepts+1
 ;     print, 'accept'
   endif else begin
      prob = exp(-0.5*(newchisq-oldchisq))
      if randomaccept(ii) lt prob then begin ;accept
         naccepts = naccepts+1
;         print, 'accept with prob', prob
      endif else begin    ;reject
         pcurrent = pold    ; do not accept the new parameters
         newchisq = oldchisq 
         nrejects = nrejects+1
 ;        print, 'reject with prob', 1.-prob
      endelse
   endelse
   ;stop
   parr(*,ii+1) = pcurrent
   chisqarr(ii+1) = newchisq
   pold = pcurrent
   oldchisq = newchisq
endfor
setplot, 14
window,0,xsize=1500,ysize=800
!p.multi = [0,4,4]
!p.charsize=2

print, 'n accepts =', naccepts
print, 'n rejects =',nrejects
print, 'acceptance ratio', float(naccepts)/float(nrejects+naccepts)



medchi = median(chisqarr)
minstep = where(chisqarr le medchi)
minstep = minstep(0)
if minstep eq 0. then minstep=1.

plot,chisqarr,xtitle='step number', ytitle='chisq',yrange=[0,3*medchi]

print, 'Trimming position is', minstep,' at chisq ', chisqarr(minstep)

chisqarr_trimmed = chisqarr
remainedind = indgen(nstep,/long)
removeind = indgen(minstep,/long)
remove, removeind,chisqarr_trimmed,remainedind
parr_trimmed = parr[*,[remainedind]]
parr_sort = 0.*parr_trimmed

nelement = n_elements(remainedind)
sigmapos = long([0.16,0.5,0.84]*nelement)
returnvalues = fltarr(3,nparam)

print,'After trimming, there are ',nelement, ' elements left.'

;Calculating outputs and plotting
for jj=0,nparam-1 do begin
   if fixparam(jj) eq 0. then goto, skip
   sortind = bsort(parr_trimmed(jj,*),asort)
   parr_Sort(jj,*) = asort
   plot,parr(jj,*),xtitle='step number', ytitle=paraname(jj),yrange=[min(parr(jj,*))-.5,max(parr(jj,*))+.5]
   ;plothist,parr_sort(jj,*),xtitle = paraname(jj),/autobin
   bin = (max(parr_sort(jj,*))-min(parr_sort(jj,*)))/100.
   if min(parr_sort(jj,*)) ne max(parr_sort(jj,*)) then begin
      plothist,parr_sort(jj,*),xhist,yhist,bin=bin,/noplot,peak=1.
      yhist= yhist/tsum(xhist,yhist)
      plot,xhist,yhist,ytitle='posterior PDF',xtitle=paraname(jj)
   endif
   print,'1 sigma range for ',paraname(jj),' is:'
   print,parr_sort(jj,[sigmapos])
   returnvalues(*,jj) =  parr_sort(jj,[sigmapos])
   
   if strmatch(name,'cswa28*') then begin
      openw,1,'cswa28result.dat',/append
      printf,1,'1 sigma range for ',paraname(jj),' is:'
      printf,1,parr_sort(jj,[sigmapos])
      close,1
   endif
   skip: continue
endfor
;stop
end
