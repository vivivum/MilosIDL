;+
; NAME:
;	CREATE_NC
;
; AUTHOR:
;	D. Orozco Suarez
;						National Astronomical Observatory of Japan,�
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE:
;	Calculate the shift and strength of the anomalous zeeman effect in ls coupling
;		and store the results "C_N" in a common block called "QUANTIC"
;	The C_N variable is a vector structure containing different lines quantic numbers
;	Example: in the case of helium the dimension is three (three lines)
;		Each dimension contains the shifts and the strengths
;		In order to create a structure more easily, the dimension of
;			C_N is set up to the bigger one
;		The number of pi and sigma components is set in the structure and
;			does not coincide with the real dimension
;			C_N(line).N_PI, C_N(line).N_SIG    (number of components)
;			C_N(line).NUP(n_pi), C_N(line).NUB(n_sig), C_N(line).NUR(n_sig)
;				shifts of principal,blue and red components
;			C_N(line).WEP(n_pi), C_N(line).WEB(n_sig), C_N(line).WER(n_sig)
;				respective strengths
;		If more than one line, the last argument is C_N(line).FO
;			the relative strenth of each line due to the log(gf)
;
; CATEGORY:
;	Splitting
;
; CALLING SEQUENCE:
;			CREATE_N_C, DATA
;
; INPUTS:
;           DATA: A scalar vector containing: the number of lines (first element);
;				the S,L, and J quantum numbers of the lower and upper level of the
;				line transition of interest; the relative strength of each line
;
; OUTPUTS:
;			NONE
;
; COMMON BLOCKS:
;			QUANTIC,C_N: store the results in this common block
;
; EXAMPLE:
;			IDL> slo1 = 2 & llo1 = 1 & jlo1 = 2   ;lower level Fe I 630.15 nm line
;			IDL> sup1 = 2 & lup1 = 2 & jup1 = 2   ;upper level
;			IDL> slo2 = 2 & llo2 = 1 & jlo2 = 1   ;lower level Fe I 630.25 nm line
;			IDL> sup2 = 2 & lup2 = 2 & jup2 = 0   ;upper level
;			IDL> rl1 = 1                          ;relative strength line 1
;			IDL> rl2 = 0.342466                   ;relative strength line 2
;			* FOR A SINGLE LINE *
;			IDL> DATA = [1.,slo1,llo1,jlo1,sup1,lup1,jup1]
;			IDL> CREATE_NC,DATA
;			* FOR TWO LINES *
;			IDL> DATA = [2.,slo1,llo1,jlo1,sup1,lup1,jup1,rl1,$
;						slo2,llo2,jlo2,sup2,lup2,jup2,rl2]
;			IDL> CREATE_NC,DATA
;
; NOTES:
;			*IMPORTANT*
;			MILOS DOES NOT WORK IF NOT INITIALIZED WITH CREATE_NC. SEE INIT_MILOS.
;
; ROUTINES CALLING COMMON BLOCK:
;           MIL_SINRF, ME_DER, MILOS
;
; CALLED ROUTINES:
;           QUANTEN, TYPE
;
; MODIFICATION HISTORY:
; 	First version created, D Orozco Suarez (DOS), 2004
;-

PRO CREATE_NC,DAT,quiet=quiet,not_normalize=not_normalize

COMMON QUANTIC,C_N

prt=keyword_set(QUIET)

a=type(c_n)
if a then c_n=0

LINES=DAT(0)     ;1,2,3 number of lines

SELECT=INDGEN(LINES)*7+1
Sl=DAT(SELECT)
Ll=DAT(SELECT+1)
Jl=DAT(SELECT+2)
Su=DAT(SELECT+3)
Lu=DAT(SELECT+4)
Ju=DAT(SELECT+5)
If Lines gt 1 then Fos=Dat(SELECT+6)

NP=0
FOR J=0,LINES-1 DO BEGIN
NPP=2.*MIN([jl(J),ju(J)])+1
IF NPP GT NP THEN NP=NPP
ENDFOR
NS=MAX(jl+ju)

NNP=FLTARR(NP)
NNS=FLTARR(NS)

;CREATE THE STRUCTURE
C_N=REPLICATE({N_PI:0.,N_SIG:0.,NUB:NNS,NUP:NNP,NUR:NNS,$
               WEB:NNS,WEP:NNP,WER:NNS,GL:0.,GU:0.,GEFF:0.,Fo:1d0},LINES)

;ASSIGN THE VARIABLES TO THE STRUCTURE
FOR I=0,LINES-1 DO BEGIN
    N_PI=2.*MIN([jl(I),ju(I)])+1
    N_SIG=jl(I)+ju(I)
    ;QUANTEN COMPUTES THE CUANTEM NUMBERS
    QUANTEN,SL(I),SU(I),LL(I),LU(I),JL(I),JU(I),$
          NUB,NUP,NUR,WEB,WEP,WER,GLO,GUP,GEF,not_normalize=not_normalize
    C_N(I).N_PI=N_PI
    C_N(I).N_SIG=N_SIG
    C_N(I).NUB(0:N_SIG-1)=NUB
    C_N(I).NUP(0:N_PI-1)=NUP
    C_N(I).NUR(0:N_SIG-1)=NUR
    C_N(I).WEB(0:N_SIG-1)=WEB
    C_N(I).WEP(0:N_PI-1)=WEP
    C_N(I).WER(0:N_SIG-1)=WER
    C_N(I).GL=GLO
    C_N(I).GU=GUP
    C_N(I).GEFF=GEF
    If Lines gt 1 then C_N(I).FO=Fos(I)

ENDFOR

IF not(prt) THEN begin
    PRINT,'------------    QUANTUM NUMBERS   --------------'
FOR I=0,LINES-1 DO BEGIN
    PRINT,'------------------------------------------------'
    PRINT,'Line:                         ',I
    PRINT,'Number of Pi components       ',C_N(I).N_PI
    PRINT,'Number of Sigma components    ',C_N(I).N_SIG
    PRINT,'Lower level lande factor      ',C_N(I).GL
    PRINT,'Upper level lande factor      ',C_N(I).GU
    PRINT,'Effective lande factor        ',C_N(I).GEFF
    PRINT,'Shifts principal component    ',C_N(I).NUP
    PRINT,'Shifts blue component         ',C_N(I).NUB
    PRINT,'Shifts red component          ',C_N(I).NUR
    PRINT,'Strength principal component  ',C_N(I).WEP
    PRINT,'Strength blue component       ',C_N(I).WEB
    PRINT,'Strength red component        ',C_N(I).WER
    PRINT,'Relative line strength        ',C_N(I).Fo
ENDFOR
    PRINT,'------------ END QUANTUM NUMBERS --------------'
endif

END
