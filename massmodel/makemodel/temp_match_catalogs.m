function temp_match_catalogs

f=load('model_colors_macs11495.txt');
g=load('model_colors_macs11495_f814w.txt');

f(:,7)=f(:,7)/10;

for i=1:length(f(:,1))
   for j=1: length(g(:,1))
    
       if f(i,2)==g(j,2) && f(i,3)==g(j,3) 
          
           f(i,7)=g(j,7);
           
       end
       
   end
end

save 'model_colors11495new2.txt' f -ascii