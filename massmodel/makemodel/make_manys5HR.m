function make_manys5HR

%s=dir('/home/wise-tvx2/data3/adiz/tmp/BestClustersSim/100Gals/HR');
%conto=0;
for q=130%:15:145
    for s=4:8:20
%        d= sprintf('q%s');
%        display(d)
if (s<10)
        strn.name(1)='q';
        strn.name(2:4)=num2str(q);
        strn.name(5)='D';
        strn.name(6)='M';
        strn.name(7)=num2str(s);
        strn.name
        mkdir(strn.name)
        cd(strn.name)
        if (s==4)
     copyfile('../qtemp/*.m','.');
     copyfile('../qtemp/*.t*','.');
     %copyfile('../../q1DM4/*.c*','.');
     go_try_gen((q/100),s)
     cd ../
        else
           older.name(1)='.';
           older.name(2)='.';
           older.name(3)='/';
         older.name(4)='q';
        older.name(5:7)=num2str(q);
        older.name(8)='D';
        older.name(9)='M';
        older.name(10)=num2str(s-8); 
        older.name(11)='/';
        older.name(12)='*';
        older.name(13)='.';
        older.name(14)='m';
        older.name(15)='*';
      copyfile(older.name,'.');
      older.name(14)='t';
      copyfile(older.name,'.');
     %copyfile('../../q1DM4/*.c*','.');
     go_try_gen2((q/100),s)
     cd ../
        end
end
if (s>12)
        strn.nameo(1)='q';
        strn.nameo(2:4)=num2str(q);
        strn.nameo(5)='D';
        strn.nameo(6)='M';
        strn.nameo(7:8)=num2str(s);
        strn.nameo
        mkdir(strn.nameo)
        cd(strn.nameo)
           older.nameo(1)='.';
           older.nameo(2)='.';
           older.nameo(3)='/';
         older.nameo(4)='q';
        older.nameo(5:7)=num2str(q);
        older.nameo(8)='D';
        older.nameo(9)='M';
        older.nameo(10:11)=num2str(s-8); 
        older.nameo(12)='/';
        older.nameo(13)='*';
        older.nameo(14)='.';
        older.nameo(15)='m';
        older.nameo(16)='*';
      copyfile(older.nameo,'.');
      older.nameo(15)='t';
      copyfile(older.nameo,'.');
     %copyfile('../../q1DM4/*.c*','.');
     go_try_gen2((q/100),s)
     cd ../
end
        
if (s==12)
        strn.namea(1)='q';
        strn.namea(2:4)=num2str(q);
        strn.namea(5)='D';
        strn.namea(6)='M';
        strn.namea(7:8)=num2str(s);
        strn.namea
        mkdir(strn.namea)
        cd(strn.namea)
           older.namea(1)='.';
           older.namea(2)='.';
           older.namea(3)='/';
         older.namea(4)='q';
        older.namea(5:7)=num2str(q);
        older.namea(8)='D';
        older.namea(9)='M';
        older.namea(10)=num2str(s-8); 
        older.namea(11)='/';
        older.namea(12)='*';
        older.namea(13)='.';
        older.namea(14)='m';
        older.namea(15)='*';
      copyfile(older.namea,'.');
      older.namea(14)='t';
      copyfile(older.namea,'.');
     %copyfile('../../q1DM4/*.c*','.');
     go_try_gen2((q/100),s)
     cd ../

        end
    end
end
