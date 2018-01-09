pro aoverb,a,upa,lowa,b,berr
;Treated as two sided gaussian
a = float(a)
upa=float(upa)
lowa=float(lowa)
b=float(b)
berr=float(berr)

f=a/b
upf=sqrt(f^2*(((upa-a)/a)^2+(berr/b)^2))
lowf=sqrt(f^2*(((a-lowa)/a)^2+(berr/b)^2))

upf = f+upf
lowf= f-lowf

print, f,upf,lowf
;stop
end
