function agnese_matchedByEye

f=load('regionsFixedbyhandNicha.txt');
g=load('temp_cat_all.txt');
%g=load('temp_cat_all_withAbmag_z0465.txt');

%h=load('temp_cat_allACSIR.txt');
counti=0;
for i=1:length(f(:,1))
    flagi(i)=0;
    for j=1:length(g(:,1))
       if sqrt( (f(i,1)-g(j,2))^2+( f(i,2)-g(j,3))^2)<3 
       display('found exact')
        counti=counti+1;
       new_file(counti,:)=g(j,:);
        flagi(i)=1;
%        elseif sqrt( (f(i,1)-h(j,2))^2+( f(i,2)-h(j,3))^2)<10 
%        display('found exact')
%         counti=counti+1;
%        new_file(counti,:)=h(j,:);
%         flagi(i)=1;
           
       end
       
    end
end


counti
display('objects found') 
a=find(flagi==0);
f(a,1:2)

% %now looking in IR only:
% for p=1:length(a)
%     i=a(p);
%     %flagi(i)=0;
%     for j=1:length(h(:,1))
%        if sqrt( (f(i,1)-h(j,2))^2+( f(i,2)-h(j,3))^2)<10 
%        display('found exact')
%         counti=counti+1;
%        new_file(counti,:)=h(j,:);
%         flagi(i)=1;
%            
%        end
%        
%     end
% end
% display('objects found') 
% a=find(flagi==0);
% f(a,1:2)
% length(a)
save model_colors_templim24updated.txt new_file -ascii