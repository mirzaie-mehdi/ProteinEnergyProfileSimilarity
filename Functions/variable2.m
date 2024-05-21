function [ v ] = variable(i)
switch (i)
    case 1
        syms FF;
        v=[FF];
    case 2
        syms LL;
        v=[LL];
    case 3
        syms II ;
        v=[II];
    case 4
        syms VV;
        v=[VV];
    case 5
        syms WW;
        v=[WW];
    case 6
        syms YY;
        v=[YY];
    case 7
        syms MM;
        v=[MM];
    case 8
        syms CC;
        v=[CC];
    case 9
        syms HH;
        v=[HH];
    case 10
        syms TT;
        v=[TT];
    case 11
        syms RR ;
        v=[RR];
    case 12
        syms AA;
        v=[AA];
    case 13
        syms NN;
        v=[NN];
    case 14
        syms QQ;
        v=[QQ];
    case 15
        syms PP;
        v=[PP];
    case 16
        zss=sym('SS');
        v=[zss];
    case 17
        syms DD;
        v=[DD];
    case 18
        syms GG;
        v=[GG];
    case 19
        syms KK;
        v=[KK];
    case 20
        syms EE;
        v=[EE];
end
