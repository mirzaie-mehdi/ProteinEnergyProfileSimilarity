 clear;clc
 tic
 path='../Data/Dunbrack/';
 % cpe=readtable([path 'cpe.csv']);
 % spe=readtable([path 'spe.csv']);
 % zi=cell(1,20);
 % x = 'FLIVWYMCHTRANQPSDGKE';
 % for NP=1:size(cpe,1)
 %     disp(num2str(NP))
 %     % ------- SPE210 to matrix -------------
 %     ms=zeros(20,20); 
 %     idx =tril(true(size(ms)));
 %     ms(idx)=table2array(spe(NP,3:212));
 %     ms=ms';
 %     ms(idx)=table2array(spe(NP,3:212));
 %     % ------- CPE210 to matrix -------------
 %     mc=zeros(20,20); 
 %     idx =tril(true(size(mc)));
 %     mc(idx)=table2array(cpe(NP,3:212));
 %     mc=mc';
 %     mc(idx)=table2array(cpe(NP,3:212));
 %     % ----------------------------------
 %     for i=1:20
 %        syms E;
 %        ei=sum(ms(i,:));
 %        for j=1:20
 %            E(j)=[x(i) x(j)];
 %            E(j)=mc(i,j)*E(j);
 %        end
 %        zi{NP,i}= (ei - sum(E))^2;
 %     end     
 % end
 % toc
 % %save([path '/ziReweight.mat'], 'zi')
 zzi=cell(1,20);
 for i=1:20
     disp(num2str(i))
     for j=2:size(zi,1)
         zi{1,i}=zi{1,i}+zi{j,i};
     end
     zzi{1,i}=zi{1,i};
 end
save([path '/zziReweight.mat'],'zzi')