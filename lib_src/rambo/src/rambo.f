
      module ol_ramboX
      contains

!------------------------------------------------------
      SUBROUTINE RAMBO(N,ET,XM,P,WT)
!------------------------------------------------------
!
!                       RAMBO
!
!    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
!
!    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
!    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
!    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
!    (MODIFIED BY R. PITTAU)
!
!                INPUT                 OUTPUT
!
!                N, ET, XM             P, WT
!
!    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
!    ET = TOTAL CENTRE-OF-MASS ENERGY
!    XM = PARTICLE MASSES ( DIM=100 )
!    P  = PARTICLE MOMENTA ( DIM=(4,100) )
!    DJ = 1/(WEIGHT OF THE EVENT)
!
!------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XM(100),P(4,100),Q(4,100),Z(100),R(4),
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(5)
      SAVE ACC,ITMAX,IBEGIN,IWARN,Z,TWOPI,PO2LOG
      DATA ACC/1.D-14/,ITMAX/10/,IBEGIN/0/,IWARN/5*0/
!
! INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      lflag= 0
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.*DATAN(1.D0)
      PO2LOG=LOG(TWOPI/4.)
      Z(2)=PO2LOG
      DO 101 K=3,100
  101 Z(K)=Z(K-1)+PO2LOG-2.*LOG(DFLOAT(K-2))
      DO 102 K=3,100
  102 Z(K)=(Z(K)-LOG(DFLOAT(K-1)))
!
! CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      STOP
!
! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+ABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP

  201 CONTINUE
      if (lflag.eq.1) then
        w0= exp((2.*N-4.)*LOG(ET)+Z(N))
        do j= 1,N
          v(j)= sqrt(p(1,j)**2+p(2,j)**2+p(3,j)**2)
        enddo

        a1= 0.d0
        a3= 0.d0
        a2= 1.d0
        do j= 1,N
          a1= a1+v(j)/ET
          a2= a2*v(j)/p(4,j)
          a3= a3+v(j)*v(j)/p(4,j)/ET
        enddo
        wm= a1**(2*N-3)*a2/a3
        dj= 1.d0/w0/wm
        return
      endif
!
! THE PARAMETER VALUES ARE NOW ACCEPTED
!
! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE

      DO 202 I=1,N
      call rans(RAN1)
      call rans(RAN2)
      call rans(RAN3)
      call rans(RAN4)
      C=2.*RAN1-1.
      S=SQRT(1.-C*C)
      F=TWOPI*RAN2
      Q(4,I)=-LOG(RAN3*RAN4)
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*COS(F)
  202 Q(1,I)=Q(4,I)*S*SIN(F)
!
! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=SQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS
!
! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
!
! CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*LOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
!
! RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
!
! MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=SQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=SQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(ABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
!
! CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*LOG(X)+LOG(WT2/WT3*ET)
!
! RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
!
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D15.6)
      END SUBROUTINE RAMBO

!-- Author :    F. James, modified by Mike Seymour
!-----------------------------------------------------------------------
      FUNCTION RANGEN(I)
!-----------------------------------------------------------------------
!     MAIN RANDOM NUMBER GENERATOR
!     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RANGEN,RANSET,RANGET
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12345,67890/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANGEN=DBLE(IZ)*4.656613001013252D-10
!--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
!-----------------------------------------------------------------------
      ENTRY RANSET(JSEED)
!-----------------------------------------------------------------------
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) then
         write(6,*) 'Jseeds=0, wrong settings for RANSET'
         stop
      endif
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
!-----------------------------------------------------------------------
      ENTRY RANGET(JSEED)
!-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END FUNCTION RANGEN

      subroutine rans(ran)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                c
!     Random number generator                                    c
!                                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision ran
      ran= rangen(1)
      end subroutine rans

      end module ol_ramboX
