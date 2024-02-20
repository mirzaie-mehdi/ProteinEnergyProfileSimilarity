 clear;clc
 tic
 path='../../Dunbrack/';
 load('../Data/rds/energy_dell_dunbrack.mat');
 Tdist=6;start=0.75;binw=0.5;seqsep =1;
 zi=cell(1,20);
 esk=cell(1,2);
 count_pair_aa=zeros(20,20);count_aa=zeros(20,1);count_pair=zeros(20,20);
 x = 'FLIVWYMCHTRANQPSDGKE';
 Nfile=1;
 fid2=fopen([path '/list_with_chainID.txt']);
 filaname_chain=fgetl(fid2);
 while ischar(filaname_chain)
     ff=lower(filaname_chain(1:4));
     disp(num2str(Nfile))
     pdb=pdbread([path 'PDBs/' ff '.pdb']);
     if length(filaname_chain)==5
         chain=filaname_chain(end);
         for i=1:length(pdb.Sequence)
             if pdb.Sequence(i).ChainID==chain
                 seq=pdb.Sequence(i).Sequence;
                 break
             end
         end
     else 
         chain=pdb.Sequence(1).ChainID;
         seq=pdb.Sequence.Sequence;
     end
     nrow=size(pdb.Model(1).Atom,2);
     coordinate=[];atnumber=[];amnumber=[];
     for i=1:nrow
         row=pdb.Model(1).Atom(i);
         if isempty(row.altLoc) || strcmp(row.altLoc,'A')
            if ~strcmp(row.AtomName(1),'H')
               if strcmp(chain,row.chainID)
                  if ~strcmp(row.AtomName,'H')&& ~strcmp(row.AtomName,'OXT')
                     coordinate=cat(1,coordinate,[row.X row.Y row.Z]);
                     atnumber=cat(1,atnumber,row.AtomSerNo);
                     amnumber=cat(1,amnumber,row.resSeq);
                  end           
               end
            end
         end
     end
     if min(amnumber)<=0
        amnumber=amnumber+abs(min(amnumber))+1; 
     end
     nr=size(coordinate,1);
     atomname=cell(nr,1);amacidname=cell(nr,1);j=1;
     for i=1:nrow
         row=pdb.Model(1).Atom(i);
         if isempty(row.altLoc)|| strcmp(row.altLoc,'A')
            if ~strcmp(row.AtomName(1),'H')
               if strcmp(chain,row.chainID)
                  if ~strcmp(row.AtomName,'H')&& ~strcmp(row.AtomName,'OXT')
                     atomname{j,1}=row.AtomName;
                     amacidname{j,1}=row.resName;j=j+1;
                  end           
               end
            end
         end
     end
     natom=numel(amnumber);
     neighbour=cell(1,natom);
     T=delaunayn(coordinate);
     try  
        tn=size(T,1);
        for k=1:tn
            for i=1:4
                X=neighbour{1,T(k,i)};
                for j=i+1:4
                    c=abs(amnumber(T(k,i))-amnumber(T(k,j)));
                    if (c>=1)
                       if isempty(find(X==T(k,j),1))
                          if ~(((strcmp(atomname{T(k,i),1},'N'))&&(strcmp(atomname{T(k,j),1},'C'))&&(c==1))||((strcmp(atomname{T(k,i),1},'C'))&&(strcmp(atomname{T(k,j),1},'N'))&&(c==1)))
                                 d=norm(coordinate(T(k,i),:)-coordinate(T(k,j),:));
                                 if d<=Tdist
                                    neighbour{1,T(k,i)}=(cat(2,neighbour{1,T(k,i)},T(k,j)));
                                    neighbour{1,T(k,j)}=(cat(2,neighbour{1,T(k,j)},T(k,i)));
                                 end
                           end
                       end
                    end
                end
            end
        end
     catch
           disp(['Error in delaunayn of ' filename])
     end
     nres=amnumber(end);
     am=zeros(nres,1);cnt=zeros(nres,nres);am_at=zeros(natom,1);ER=zeros(nres,nres);atnumber=zeros(natom,1);
     int=[];
     for i=1:natom
         Nei=neighbour{1,i}; nNei=numel(Nei);
         ami=amacidname{i};am_at(i,1)=amacid2num(ami);
         try
            ssi=amnumber(i);am(ssi,1)=amacid2num(ami);atnumber(i,1)=ssi;count_aa(am(ssi,1),1)=count_aa(am(ssi,1),1)+1;           
            crj=coordinate(i,:);
            ati=atomname{i,1};
             vi=atomtype2num167(ati,ami);
             for k=1:nNei
                 ssk=amnumber(Nei(k));
                 ss=(abs(ssi-ssk));
                 amk=amacidname{Nei(k)};
                 crk=coordinate(Nei(k),:);
                 atk=atomname{Nei(k)};
                 vk=atomtype2num167(atk,amk);
                 d=norm(crj-crk);
                 if d<=Tdist && ss>=seqsep && vi>0 && vk>0
                    F=fix((d-start)/binw)+1;ss=1;F=max(F,1);s1=min(vi,vk);t1=max(vi,vk);
                    if  ~isinf(energy_dell_dunbrack{s1,t1}(ss,F)) 
                        ER(ssi,ssk)=ER(ssi,ssk)+energy_dell_dunbrack{s1,t1}(ss,F);
                        ER(ssk,ssi)=ER(ssk,ssi)+energy_dell_dunbrack{s1,t1}(ss,F);
                        cnt(ssi,ssk)=cnt(ssi,ssk)+1;cnt(ssk,ssi)=cnt(ssk,ssi)+1;
                        atamik=cell(1,4);
                        atamik{1,1}=ati;atamik{1,2}=atk;atamik{1,3}=ami;atamik{1,4}=amk;
                        int=cat(1,int,atamik);
                    else
                        disp([num2str(s1) ',' num2str(t1) ',' num2str(F)])
                    end
                 end
             end
         catch        
         end 
     end   
     es=zeros(20,20);  
     %count_pair=zeros(20,20); 
     try
        for u=1:20
            a1=find(am==u);
            for s=1:20
                a2=find(am==s);
                %count_pair(u,s)=sum(sum(cnt(a1,a2)));
                es(u,s) =sum(sum(ER(a1,a2)));
            end
        end
     catch
     end
     esk{Nfile,1}=es;
     esk{Nfile,2}=nres;
     ni= zeros(1,20);
     for i=1:20
          pos = strfind(seq,x(i));
          ni(i) = size(pos,2);
     end 
     %ni = ni/length(seq);
     for i=1:20
        syms E;
        ei=sum(es(i,:));
        for j=1:20
            E(j)=[x(i) x(j)];
            E(j)=(ni(j)/length(seq))*E(j);
        end
        zi{Nfile,i}= (ei - (ni(i)*sum(E)))^2;
     end
     %%%%%%%%%%%%%
     filaname_chain=fgetl(fid2); 
     Nfile=Nfile+1;      
 end
 toc
 save([path '/zi.mat'], 'zi')
 zzi=cell(1,20);
 for i=1:20
     for j=2:size(zi,1)
         zi{1,i}=zi{1,i}+zi{j,i};
     end
     zzi{1,i}=zi{1,i};
 end
save([path '/zzi.mat'],'zzi')
Pij = zeros(20,20);
for i=1:20
    disp(num2str(i))
    v = variable(i);
    grd = gradient(zzi{1,i},v);
    slv = solve(grd);
    slv_cel= struct2cell(slv);
    slv_cel2=slv_cel([5 10 8 18 19 20 11 2 7 17 15 1 12 14 13 16 3 6 9 4]);
    for j=1:20
        Pij(i,j) = vpa(slv_cel2{j},10);
    end
end
writematrix(Pij,[path '/Pij.csv'])