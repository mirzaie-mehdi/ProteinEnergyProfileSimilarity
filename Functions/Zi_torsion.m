 clc
 tic
 path='../Data/Dunbrack/';
 zi=cell(1,20);
 x = 'FLIVWYMCHTRANQPSDGKE';
 for NP=1:6384
     disp(num2str(NP))
     l=len(NP);
     for i=1:20
         ni=freqAA20(NP,i);
         ni = ni/l;
         syms E;
         E(1)=[x(i) x(i)];
         ei=E20torsion(NP,i);
         zi{NP,i}= (ei - (ni*E))^2;
     end     
 end
 toc
 save([path '/zi_torsion.mat'], 'zi')
 zzi=cell(1,20);
 for i=1:20
     for j=2:size(zi,1)
         zi{1,i}=zi{1,i}+zi{j,i};
     end
     zzi{1,i}=zi{1,i};
 end
save([path '/zzi_torsion.mat'],'zzi')