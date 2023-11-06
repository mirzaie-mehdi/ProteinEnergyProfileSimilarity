function [ v ] = variable(i)
switch (i)
    case 1
        zfi=sym('FI');
        syms FF FL FV FW FY FM FC FH FT FR FA FN FQ FP FS FD FG FK FE;
        v=[FF FL zfi FV FW FY FM FC FH FT FR FA FN FQ FP FS FD FG FK FE];
    case 2
        zls=sym('LS');
        syms LF LL LI LV LW LY LM LC LH LT LR LA LN LQ LP LD LG LK LE;
        v=[LF LL LI LV LW LY LM LC LH LT LR LA LN LQ LP zls LD LG LK LE];
    case 3
        zif=sym('IF');
        ziv=sym('IV');
        syms IL II IW IY IM IC IH IT IR IA IN IQ IP IS ID IG IK IE;
        v=[zif IL II ziv IW IY IM IC IH IT IR IA IN IQ IP IS ID IG IK IE];
    case 4
        syms VF VL VI VV VW VY VM VC VH VT VR VA VN VQ VP VS VD VG VK VE;
        v=[VF VL VI VV VW VY VM VC VH VT VR VA VN VQ VP VS VD VG VK VE];
    case 5
        syms WF WL WI WV WW WY WM WC WH WT WR WA WN WQ WP WS WD WG WK WE;
        v=[WF WL WI WV WW WY WM WC WH WT WR WA WN WQ WP WS WD WG WK WE];
    case 6
        syms YF YL YI YV YW YY YM YC YH YT YR YA YN YQ YP YS YD YG YK YE;
        v=[YF YL YI YV YW YY YM YC YH YT YR YA YN YQ YP YS YD YG YK YE];
    case 7
        syms MF ML MI MV MW MY MM MC MH MT MR MA MN MQ MP MS MD MG MK ME;
        v=[MF ML MI MV MW MY MM MC MH MT MR MA MN MQ MP MS MD MG MK ME];
    case 8
        zcv=sym('CV');
        syms CF CL CI CW CY CM CC CH CT CR CA CN CQ CP CS CD CG CK CE;
        v=[CF CL CI zcv CW CY CM CC CH CT CR CA CN CQ CP CS CD CG CK CE];
    case 9
        syms HF HL HI HV HW HY HM HC HH HT HR HA HN HQ HP HS HD HG HK HE;
        v=[HF HL HI HV HW HY HM HC HH HT HR HA HN HQ HP HS HD HG HK HE];
    case 10
        ztf=sym('TF');
        syms TL TI TV TW TY TM TC TH TT TR TA TN TQ TP TS TD TG TK TE;
        v=[ztf TL TI TV TW TY TM TC TH TT TR TA TN TQ TP TS TD TG TK TE];
    case 11
        syms RF RL RI RV RW RY RM RC RH RT RR RA RN RQ RP RS RD RG RK RE;
        v=[RF RL RI RV RW RY RM RC RH RT RR RA RN RQ RP RS RD RG RK RE];
    case 12
        zar=sym('AR');
        syms AF AL AI AV AW AY AM AC AH AT AA AN AQ AP AS AD AG AK AE;
        v=[AF AL AI AV AW AY AM AC AH AT zar AA AN AQ AP AS AD AG AK AE];
    case 13
        syms NF NL NI NV NW NY NM NC NH NT NR NA NN NQ NP NS ND NG NK NE;
        v=[NF NL NI NV NW NY NM NC NH NT NR NA NN NQ NP NS ND NG NK NE];
    case 14
        syms QF QL QI QV QW QY QM QC QH QT QR QA QN QQ QP QS QD QG QK QE;
        v=[QF QL QI QV QW QY QM QC QH QT QR QA QN QQ QP QS QD QG QK QE];
    case 15
        zpy=sym('PY');
        zpe=sym('PE');
        syms PF PL PI PV PW PM PC PH PT PR PA PN PQ PP PS PD PG PK;
        v=[PF PL PI PV PW zpy PM PC PH PT PR PA PN PQ PP PS PD PG PK zpe];
    case 16
        zsf=sym('SF');
        zss=sym('SS');zsl=sym('SL');
        syms SI SV SW SY SM SC SH ST SR SA SN SQ SP SD SG SK SE;
        v=[zsf zsl SI SV SW SY SM SC SH ST SR SA SN SQ SP zss SD SG SK SE];
    case 17
        syms DF DL DI DV DW DY DM DC DH DT DR DA DN DQ DP DS DD DG DK DE;
        v=[DF DL DI DV DW DY DM DC DH DT DR DA DN DQ DP DS DD DG DK DE];
    case 18
        zgf=sym('GF');
        zga=sym('GA');
        syms GL GI GV GW GY GM GC GH GT GR GN GQ GP GS GD GG GK GE;
        v=[zgf GL GI GV GW GY GM GC GH GT GR zga GN GQ GP GS GD GG GK GE];
    case 19
        syms KF KL KI KV KW KY KM KC KH KT KR KA KN KQ KP KS KD KG KK KE;
        v=[KF KL KI KV KW KY KM KC KH KT KR KA KN KQ KP KS KD KG KK KE];
    case 20
        syms EF EL EI EV EW EY EM EC EH ET ER EA EN EQ EP ES ED EG EK EE;
        v=[EF EL EI EV EW EY EM EC EH ET ER EA EN EQ EP ES ED EG EK EE];
end
