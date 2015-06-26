;+
; NAME: 
;	QUANTEN
;
; AUTHOR: 
;	D. Orozco Suarez  	
;						National Astronomical Observatory of Japan,Ê
;						2-21-1 Osawa, Mitaka, 181-8588, JAPAN
;						d.orozco@nao.ac.jp
;
; PURPOSE: 
;	Evaluate the atomic quanten numbers for the anomalous zeeman effect in LS coupling
;
; CATEGORY: 
;	Splitting
;
; CALLING SEQUENCE: 
;			QUANTEN, S1,S2,L1,L2,J1,J2,msig1,mpi,msig2,sig1,pi,sig2,g1,g2,geff
;
; INPUTS:
;			S1,S2,L1,L2,J1,J2: S, L, and J quantum numbers of the line transition 
;				for the lower (...1) and upper (...2) level 
; OUTPUTS:
;			msig1, mpi, msig2, sig1, pi, sig2: position and strength of zeeman components
;			g1, g2: Level landŽ factors
;			geff; Effective LandŽ factor
;           
; COMMON BLOCKS:
;			NONE 
;
; CALLED ROUTINES:
;           NONE
;
; MODIFICATION HISTORY:
; 	First version created, D Orozco Su‡rez (DOS), 2004
;
;-


PRO QUANTEN,S1,S2,L1,L2,J1,J2,msig1,mpi,msig2,sig1,pi,sig2,g1,g2,geff

;1-> LOW
;2-> UP

;magnetic quanten number Mlo Mup
m1=indgen(2*j1+1)-j1
m2=indgen(2*j2+1)-j2
JJ = J2-J1
N_PI=2.*MIN([j1,j2])+1 ;Number of pi components
N_SIG=j1+j2 ;Number of sigma components

;lande factors (with 'if' because j could be cero)
if j1 ne 0 then g1=(3.*j1*(j1+1)+s1*(s1+1)-l1*(l1+1))/(2.*j1*(j1+1)) else g1=0.
;geff=(g1+g2)/2.+(g1-g2)*d/4.
if j2 ne 0 then g2=(3.*j2*(j2+1)+s2*(s2+1)-l2*(l2+1))/(2.*j2*(j2+1)) else g2=0.
;geff=(g1+g2)/2.+(g1-g2)*d/4.
geff=Double((g1+g2)/2.+(g1-g2)*(j1*(j1+1.)-j2*(j2+1.))/4.)


;BLUE COMPONENT => Mlo-Mup = +1
;RED COMPONENT => Mlo-Mup = -1
;CENTRAL COMPONENT => Mlo-Mup = 0
pi=FLTARR(N_PI)
sig1=FLTARR(N_SIG)
sig2=FLTARR(N_SIG)
mpi=FLTARR(N_PI)
msig1=FLTARR(N_SIG)
msig2=FLTARR(N_SIG)

;COUNTERS FOR THE COMPONENTS
ipi=0
isig1=0
isig2=0

for j=0,2*j1 do begin
    for i=0,2*j2 do begin
        IM=M2(i)-M1(j)
        case IM of
            0: begin            ;M -> M  ;CENTRAL COMPONENT
                case JJ of
                    -1: begin   ;  j -> j-1
                        pi(ipi)=j1^2.-m1(j)^2.
                        mpi(ipi)=g1*m1(j)-g2*m2(i)
                        ipi=ipi+1
                    end
                    0: begin    ;  j -> j
                        pi(ipi)=m1(j)^2.
                        mpi(ipi)=g1*m1(j)-g2*m2(i)
                        ipi=ipi+1
                    end
                    1: begin    ;  j -> j+1
                        pi(ipi)=(j1+1.)^2.-m1(j)^2.
                        mpi(ipi)=g1*m1(j)-g2*m2(i)
                        ipi=ipi+1
                    end
                    else: print,''
                endcase
            end
            1: begin            ;M -> M+1  ;BLUE COMPONENT
                case JJ of
                    -1: begin   ;  j -> j-1
                        sig1(isig1)=(j1-m1(j))*(j1-m1(j)-1)/4.
                        msig1(isig1)=g1*m1(j)-g2*m2(i)
                        isig1=isig1+1
                    end
                    0: begin    ;  j -> j
                        sig1(isig1)=(j1-m1(j))*(j1+m1(j)+1)/4.
                        msig1(isig1)=g1*m1(j)-g2*m2(i)
                        isig1=isig1+1
                    end
                    1: begin    ;  j -> j+1
                        sig1(isig1)=(j1+m1(j)+1)*(j1+m1(j)+2)/4.
                        msig1(isig1)=g1*m1(j)-g2*m2(i)
                        isig1=isig1+1
                    end
                    else: print,''
                endcase
            end
            -1: begin           ;M -> M-1   ;RED COMPONENT
                case JJ of
                    -1: begin   ;  j -> j-1
                        sig2(isig2)=(j1+m1(j))*(j1+m1(j)-1)/4.
                        msig2(isig2)=g1*m1(j)-g2*m2(i)
                        isig2=isig2+1
                    end
                    0: begin    ;  j -> j
                        sig2(isig2)=(j1+m1(j))*(j1-m1(j)+1)/4.
                        msig2(isig2)=g1*m1(j)-g2*m2(i)
                        isig2=isig2+1
                    end
                    1: begin    ;  j -> j+1
                        sig2(isig2)=(j1-m1(j)+1)*(j1-m1(j)+2)/4.
                        msig2(isig2)=g1*m1(j)-g2*m2(i)
                        isig2=isig2+1
                    end
                    else: print,''
                endcase
            end
            else: goto,salto ;NOT ALLOWED TRANSITION
        endcase
        salto:
    endfor
endfor

;normalization OF EACH COMPONENT

pi=pi/total(pi)
sig1=sig1/total(sig1)
sig2=sig2/total(sig2)

end


