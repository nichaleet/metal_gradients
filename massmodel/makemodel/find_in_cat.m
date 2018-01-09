function find_in_cat(x,y)

f=load('temp_cat_all.txt');
for i=1:length(f(:,1))
   if abs(f(i,2)-x)<8 && abs(f(i,3)-y)<8
      i 
      f(i,:)
   end
    
end