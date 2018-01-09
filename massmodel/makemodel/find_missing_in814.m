function find_missing_in814(x,y)
f=load('814byAdi.cat');
lines=length(f(:,1));
for i=1:lines
    if (x<f(i,2)+15 && x>f(i,2)-15 && y<f(i,3)+15 && y>f(i,3)-15)
        f(i,21:23)
        i
    end
end
%save 'tempmiss.txt' g -ascii