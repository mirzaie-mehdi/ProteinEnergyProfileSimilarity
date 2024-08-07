function [v]=atomtype2num167(atname,amname)
switch (amname)
       case 'GLY'   %*****************GLY***********GLY***********GLY****************
        switch(atname)
            case 'N'
                v=1;
            case 'CA'            %*****************backbone atoms
                v=2;
            case 'CB'
                v=2;
            case 'C'
                v=3;
            case 'O'
                v=4;
            otherwise
                v=0;
        end
    case 'ALA'                    %********ALA*****************ALA************ALA*********************
        switch (atname)
            case 'N'
                v=5;
            case 'CA'       %edited
                v=6;
            case 'C'
                v=7;
            case 'O'
                v=8;
            case 'CB'
                v=9;
            otherwise
                v=0;
        end
    case 'VAL'      %****************VAL************VAL***********VAL*********
        switch(atname)
            case 'N'
                v=10;
             case 'CA'            %*****************backbone atoms
                v=11;
            case 'C'
                v=12;
            case 'O'
                v=13;
            case 'CB'
                v=14;
            case{'CG1'}
                v=15;
            case 'CG2'
                v=16;
            otherwise
                v=0;
        end    
    case 'LEU'    %******LEU***********LEU***********LEU***************
        switch(atname)
            case 'N'
                v=17;
            case 'CA'            %*****************backbone atoms
                v=18;
            case 'C'
                v=19;
            case 'O'
                v=20;
            case 'CB'
                v=21;
            case'CG'
                v=22;
            case {'CD1'}
                v=23;
            case 'CD2'
                v=24;
            otherwise
                v=0;
        end
     case 'ILE'     %***********ILE***********ILE************ILE**************
        switch(atname) 
            case 'N'
                v=25;
            case 'CA'            %*****************backbone atoms
                v=26;
            case 'C'
                v=27;
            case 'O'
                v=28;
            case 'CB'
                v=29;
            case 'CG1' 
                v=30;
            case 'CG2'
                v=31;
            case {'CD1'}
                v=32;
            otherwise
                v=0;
        end
     case 'MET'     %**************MET******************MET*************MET*************
        switch(atname)
             case 'N'
                v=33; 
            case 'CA'            %*****************backbone atoms
                v=34;                
            case 'C'
                v=35;
            case 'O'
                v=36;
            case 'CB'
                v=37;
            case 'CG'
                v=38;
            case 'SD' 
                v=39;
            case 'CE'
                v=40;
            otherwise
                v=0;
        end
       case 'PRO'   %***********PRO************PRO**************PRO***************
        switch (atname)
            case 'N'
                v=41;    
            case 'CA'            %*****************backbone atoms
                v=42;
            case 'C'
                v=43;
            case 'O'
                v=44;
            case {'CB'}
                v=45;
             case {'CG'}
                v=46;
            case{'CD'}
                v=47;
            otherwise
                v=0;
        end
    case 'HIS'   %***************HIS***********HIS**********HIS*************
        switch(atname) 
            case 'N'
                v=48;
            case 'CA'            
                v=49;
            case 'C'
                v=50;
            case 'O'
                v=51;
            case 'CB'
                v=52;
            case 'CG'
                v=53;
            case{'ND1'}%??????????????
                v=54;
            case 'CD2'
                v=55;
            case{'CE1'}
                v=56;
             case 'NE2'%?????????????
                v=57;
              otherwise
                v=0;
        end
   case 'PHE'    %************PHE************PHE***********PHE**************
        switch(atname)
            case 'N'
                v=58;
            case 'CA'            %*****************backbone atoms
                v=59;
            case 'C'
                v=60;
            case 'O'
                v=61;
            case 'CB'
                v=62;
            case 'CG'
                v=63;
            case'CD1'
              v=64;
            case 'CD2'
                v=65;
            case 'CE1'
                v=66;
            case 'CE2'
                v=67;
            case 'CZ'
                v=68;
            otherwise
                v=0;
        end
     case 'TYR'            %**************TYR*************TYR***********TYR************
        switch(atname)
            case 'N'
                v=69;
            case 'CA'            %*****************backbone atoms
                v=70;
            case 'C'
                v=71;
            case 'O'
                v=72;
            case 'CB'
                v=73;
            case 'CG'
                v=74;
            case{'CD1'}
                v=75;
            case 'CD2'
                v=76;
            case 'CE1'
                v=77;
            case 'CE2'
                v=78;
            case {'CZ'}
                 v=79;
            case 'OH'   
                v=80;
            otherwise
                v=0;
        end
    case 'TRP'         %*********TRP**************TRP*****************TRP*************   
        switch(atname)
            case 'N'
                v=81;
            case 'CA'            %*****************backbone atoms
                v=82;
            case 'C'
                v=83;
            case 'O'
                v=84;
            case 'CB'
                v=85;
            case'CG'
                v=86;
            case 'CD1'
                v=87;
            case 'NE1'
                v=88;                
            case 'CD2'%???????????
                v=89;
            case {'CE2'}
                v=90;
            case 'CE3'
                v=91;                
            case 'CZ3'
                v=92;
            case 'CZ2'
                v=93;
            case 'CH2'
                v=94;
            otherwise
                v=0;
        end
    case{'CYH'}  %**************CYS****************CYS************CYS**********
        switch(atname)
            case 'N'
                v=95;
            case 'CA'            %*****************backbone atoms
                v=96;
            case 'C'
                v=97;
            case 'O'
                v=98;
            case 'CB'
                v=99;
            case{'SG'}
                v=100;
            otherwise
                v=0;
        end    
    case{'CYS','CYSS'}  %**************CYS****************CYS************CYS**********
        switch(atname)
            case 'N'
                v=95;
            case 'CA'            %*****************backbone atoms
                v=96;
            case 'C'
                v=97;
            case 'O'
                v=98;
            case 'CB'
                v=99;
            case{'SG'}
                v=100;
            otherwise
                v=0;
        end
    case 'SER' %*****************SER************SER*************SER*************
        switch(atname)
            case 'N'
                v=101;
            case 'CA'            %*****************backbone atoms
                v=102;
            case 'C'
                v=103;
            case 'O'
                v=104;
            case 'CB'
                v=105;
            case 'OG'   
                v=106; 
            otherwise
                v=0;
        end
    case 'THR'     %****************THR*******THR**************THR********************
        switch(atname)
             case 'N'
                v=107;
            case 'CA'            %*****************backbone atoms
                v=108;
            case 'C'
                v=109;
            case 'O'
                v=110;
            case 'CB'
                v=111;
            case 'OG1' 
                v=112;
            case 'CG2'  
                v=113;
            otherwise
                v=0;
        end
    case 'ASN'        %*********ASN***********ASN*******ASN***********ASN*****************
        switch(atname)
            case 'N'
                v=114;
            case 'CA'            %*****************backbone atoms
                v=115;
            case 'C'
                v=116;
            case 'O'
                v=117;
            case 'CB'
                v=118;
            case 'CG'
                v=119;
            case 'OD1' 
                v=120;
            case 'ND2'
                v=121;
            otherwise
                v=0;
        end
   case 'GLN'        %********GLN************GLN***************GLN****************
        switch(atname)  
            case 'N'
                v=122;
            case 'CA'            %*****************backbone atoms
                v=123;
            case 'C'
                v=124;
            case 'O'
                v=125;
            case {'CB'}
                v=126;
            case 'CG'
                v=127; 
            case{'CD'}
                v=128;  
            case{'OE1'}
                v=129; 
            case{'NE2'}
                v=130;
            otherwise
                v=0;
        end
   case{'ASP'}       %**********ASP**************ASP******************ASP*********
        switch(atname)
            case 'N'
                v=131;
            case 'CA'            %*****************backbone atoms
                v=132;
            case 'C'
                v=133;
            case 'O'
                v=134;
            case 'CB'
                v=135;
            case{'CG'}
                v=136; 
            case{'OD1'}
                v=137;
            case 'OD2'
                v=138; 
            otherwise
                v=0;
        end
    case{'GLU'} %**************GLU*************GLU************GLU****************
        switch(atname)
            case 'N'
                v=139;
            case 'CA'            %*****************backbone atoms
                v=140;
            case 'C'
                v=141;
            case 'O'
                v=142;
            case {'CB'}
                v=143;
            case 'CG'
                v=144;
            case{'CD'}
                v=145;
            case{'OE1'}
                v=146;
            case 'OE2'
                v=147;
            otherwise
                v=0;
        end
    case 'LYS'   %************LYS***********LYS************LYS***************
        switch(atname)            
             case 'N'
                v=148;
            case 'CA'            %*****************backbone atoms
                v=149;
            case 'C'
                v=150;
            case 'O'
                v=151;      
            case {'CB'}
                v=152;
            case 'CG'
                v=153;
            case 'CD'
                v=154;
            case{'CE'}
                v=155;
            case 'NZ'
                v=156;
            otherwise
                v=0;
        end
  
    case 'ARG'                       %*********************ARG*********ARG***************ARG**************
        switch(atname)    
            case 'N'
                v=157;
            case 'CA'            %*****************backbone atoms
                v=158;
            case 'C'
                v=159;
            case 'O'
                v=160;
            case {'CB'}
                v=161;
            case{'CG'}
                v=162;
            case'CD'
                v=163;
            case 'NE'
                v=164;
            case 'CZ'
                v=165;
            case {'NH1'}
                v=166;
            case {'NH2'}
                v=167;
             otherwise
                v=0;
        end 
end
