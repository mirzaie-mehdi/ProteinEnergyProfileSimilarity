%   Compute the knowledge-based potential. 
%   Author: "Peyman & Mehdi"
% -------------------------------------------
 clear
 clc
 rcsb = 'https://files.rcsb.org/download/';
% path = '~/Desktop/proj/Dunbrack/PDBs/';
 Tdist=15;
 start=0.75;binw=0.5;
 Tind=fix((Tdist-start)/binw)+1;
 pair_freq=cell(167,167);
 for i=1:167
     for j=i:167
         pair_freq{i,j}=zeros(1,Tind);
     end
 end
 Nfile=1;
 fid2=fopen('list_with_chainID_rm_Olaps.txt');
 filaname_chain=fgetl(fid2);
 while ischar(filaname_chain)
     chain=filaname_chain(end);
     ff=lower(filaname_chain(1:4));
     disp([num2str(Nfile) ',' ff])
     temp_file = websave('temp.pdb', [rcsb ff '.pdb']);
     pdb=pdbread(temp_file);
%     pdb=pdbread([path ff '.pdb']);
     nrow=size(pdb.Model(1).Atom,2);
     coordinate=[];atnumber=[];amnumber=[];
     for i=1:nrow
         row=pdb.Model(1).Atom(i);
         if strcmp(chain,row.chainID)
             if ~strcmp(row.AtomName,'H')&& ~strcmp(row.AtomName,'OXT')
                 coordinate=cat(1,coordinate,[row.X row.Y row.Z]);
                 atnumber=cat(1,atnumber,row.AtomSerNo);
                 amnumber=cat(1,amnumber,row.resSeq);
             end           
         end
     end
     nr=size(coordinate,1);
     atomname=cell(nr,1);amacidname=cell(nr,1);j=1;
     for i=1:nrow
         row=pdb.Model(1).Atom(i);
         if strcmp(chain,row.chainID)
             if ~strcmp(row.AtomName,'H')&& ~strcmp(row.AtomName,'OXT')
                atomname{j,1}=row.AtomName;
                amacidname{j,1}=row.resName;j=j+1;
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
     for j=1:natom
         A=neighbour{1,j};n=numel(A);
         try
             s0=atomtype2num167(atomname{j},amacidname{j});
             crj=coordinate(j,:);
             for k=1:n
                 t0=atomtype2num167(atomname{A(k)},amacidname{A(k)});
                 s=min(s0,t0);t=max(s0,t0);
                 crk=coordinate(A(k),:);
                 ds=norm(crj-crk);
                 dd=max(fix((ds-start)/binw)+1,1);
                 if  s>0 && t>0         
                     pair_freq{s,t}(1,dd)=pair_freq{s,t}(1,dd)+1;
                 end
             end
         catch
         end
     end   
     filaname_chain=fgetl(fid2); 
     ff=lower(filaname_chain);
     Nfile=Nfile+1;
 end
save('pair_freq_rm_Olaps.mat', 'pair_freq')
% ----------------------------------------------
fr=zeros(1,29);
natom=167;
pr=cell(natom,natom);
energy_dell_dunbrack=cell(natom,natom);
M=zeros(natom,natom);
for i=1:natom
    for j=i:natom
        fre=pair_freq{i,j}+eps;
        fr=fr+fre;
        pr{i,j}=fre/sum(fre);
        M(i,j)=sum(fre);
    end
end
pfr=fr./sum(fr);
for i=1:natom
    for j=i:natom
        energy_dell_dunbrack{i,j}=0.582*log(1+M(i,j)*(1/50))+log(1+M(i,j)*(1/50)*(pr{i,j}./pfr))*(-0.582);
    end
end
save('energy_dell_dunbrack_rm_Olaps.mat', 'energy_dell_dunbrack')
% ------------------------------------------------------------------
% --------------------- Note:
% ---- Later The `energy_dell_dunbrack.mat` file was converted to the CSV file and save at `Data/csv/energy.csv`.
% ------------------------------------------------------------------
