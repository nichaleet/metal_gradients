function make_joint_catalog2;
f=load('temp_cat_all_newcoords006.txt');
g=load('A2744MembersSpecs.txt');
lines1=length(f(:,1));
lines2=length(g(:,1));
lines1
for i=1:lines1
    i
    new_file(i,1:14)=f(i,1:14);
    objflag=0;
   for j=1:lines2
    if ( abs(f(i,11)-g(j,1))<0.0005 && abs(f(i,12)-g(j,2))<0.002 ) %&& (file_data(obj,6)>=15.5) && (file_data(obj,6)<=23))
  objflag=objflag+1;
 
        new_file(i,15)=g(j,4)/(29979245800/100/1000);

    end
   end
   if objflag>1
      display('found more than 1 match for obj') 
   elseif objflag<1
       new_file(i,15)=-99;
   end
end
%    a=find(new_file(:,2)~=0);
%    new_file2=new_file(a,:);
save temp_cat_all_newcoords006wspecscol14.txt new_file -ascii;

lines1
%file_data=load('model_colors.txt');