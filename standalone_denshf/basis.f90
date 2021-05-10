 Module HObasis

    Use hfodd_types
    Use hfodd_sizes

    Implicit None

     ! Public variables broadcast to the world
    Integer(ipr) :: NXMAXX,NYMAXX,NZMAXX,NXHERM,NYHERM,NZHERM,LDBASE
    Integer(ipr), Allocatable :: NXVECT(:),NYVECT(:),NZVECT(:)
    Integer(ipr), Allocatable :: INDICE(:,:,:)
    Integer(ipr), Allocatable :: LAXOFY(:),LAXOFZ(:),LAYOFZ(:),LAYOFX(:),LAZOFX(:),LAZOFY(:)
    Integer(ipr), Allocatable :: LAXOYZ(:,:),LAYOZX(:,:),LAZOXY(:,:)
    Real(pr), Allocatable, PUBLIC :: WGTACT(:,:),PNTACT(:,:),EXPAUX(:,:)
    Real(pr), Allocatable, PUBLIC :: HERACT(:,:,:),DHRACT(:,:,:),DDHACT(:,:,:)
    Complex(pr), Allocatable, PUBLIC :: CERMTS(:,:,:),CHRMTS(:,:,:),CDHRMT(:,:,:)

    ! Private variables local to this module
    Integer(ipr), PRIVATE :: number_shells, number_HOwave, number_Cartesian
    Integer(ipr), PRIVATE :: number_Gauss, number_Gauss_max
    Integer(ipr), PRIVATE :: number_Gauss_x, number_Gauss_y, number_Gauss_z
    Integer(ipr), Allocatable, PRIVATE :: integration_grid(:)

    Real(pr) :: HBAROX=8.0_pr,HBAROY=8.0_pr,HBAROZ=8.0_pr
    Real(pr) :: scaling,precision_Hermite=1.E-14
    Real(pr), Allocatable :: HOMSCA(:)
    Real(pr), Allocatable :: HERFAC(:),FACTOR(:),oscillator_scaling(:)

 Contains

    !---------------------------------------------------------------------!
    ! This subroutine calculates harmonic oscillator wave functions and   !
    ! their first and second derivatives, as well as prepares points and  !
    ! weights for the Gauss-Hermite integration. Exponential factors are  !
    ! *NOT* included neither in the wave functions nor in the derivatives.!
    ! It provides a unique interface to routine build_HOwavefunctions,    !
    ! which does the actual job for coordinate i_Cartesian                !
    !---------------------------------------------------------------------!
    Subroutine DEFINT(NXHERM,NYHERM,NZHERM,NOBODY,HOSCAX,HOSCAY,HOSCAZ)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM,NOBODY
      Real(pr), INTENT(IN) :: HOSCAX,HOSCAY,HOSCAZ

      Integer(ipr) :: i_Cartesian

      If(NOBODY.Gt.2) Stop 'NOBODY TOO  BIG  IN DEFINT'
      If(NOBODY.Lt.1) Stop 'NOBODY TOO SMALL IN DEFINT'

      number_shells = NDOSCI
      number_HOwave = 2*NDOSCI+2
      number_Cartesian = NDKART

      Allocate(oscillator_scaling(1:number_Cartesian))
      oscillator_scaling(1:number_Cartesian) = (/ HOSCAX,HOSCAY,HOSCAZ /)

      Allocate(integration_grid(1:number_Cartesian))
      number_Gauss_x=NXHERM; integration_grid(1)=number_Gauss_x
      number_Gauss_y=NYHERM; integration_grid(2)=number_Gauss_y
      number_Gauss_z=NZHERM; integration_grid(3)=number_Gauss_z

      number_Gauss_max = Max(number_Gauss_x,number_Gauss_y)
      number_Gauss_max = Max(number_Gauss_z,number_Gauss_max)
      Allocate(FACTOR(0:number_Gauss_max))

      ! If NOBODY=2 we scale the weights and the zeros by 1/sqrt(2). This
      ! scaling corresponds to integrating the products of *FOUR* harmonic
      ! oscillator wave functions. Except for the density-dependent term of the
      ! Skyrme interaction and the Coulomb exchange term in the mean field,
      ! total energies and matrix elements of the one-body potential contain
      ! products of *FOUR* wave functions. The scaled weights and zeros will be
      ! used to integrate these quantities.
      !
      ! if NOBODY=1 no sqrt(2) scaling; in this case the weights and zeros can
      ! be used to integrate products of *TWO* harmonic oscillator wave
      !functions.
      If(NOBODY.Eq.2) scaling=Sqrt(2.0_pr)
      If(NOBODY.Eq.1) scaling=1.0_pr

      If(.Not.Allocated(WGTACT)) Allocate(WGTACT(1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(PNTACT)) Allocate(PNTACT(1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(EXPAUX)) Allocate(EXPAUX(1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(HERACT)) Allocate(HERACT(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(DHRACT)) Allocate(DHRACT(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(DDHACT)) Allocate(DDHACT(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))

      ! Define HO wave function normalization factors
      Call NORHER

      ! Defining factorials
      Call DEFFAC

      ! Build HO wave functions and their derivative on the Gauss-Hermite mesh
      Do i_Cartesian=1,number_Cartesian
         number_Gauss = integration_grid(i_Cartesian)
         Call build_HO_wavefunctions(i_Cartesian)
         Call build_HO_exponentials(i_Cartesian)
      End Do

      Deallocate(oscillator_scaling)
      Deallocate(integration_grid)
      Deallocate(FACTOR,HERFAC)

    End Subroutine DEFINT

    !---------------------------------------------------------------------!
    ! This subroutine calculates harmonic oscillator wave functions and   !
    ! their first and second derivatives, as well as prepares points and  !
    ! weights for the Gauss-Hermite integration. Exponential factors are  !
    ! *NOT* included neither in the wave functions nor in the derivatives.!
    !---------------------------------------------------------------------!
    Subroutine build_HO_wavefunctions(i_Cartesian)

      Integer(ipr), INTENT(IN) :: i_Cartesian

      Integer(ipr) :: i_Gauss,i_wave
      Real(pr) :: XZER

      Real(pr), Allocatable :: PHERMI(:),DHERMI(:)
      Real(pr), Allocatable :: XHERMI(:),WHERMI(:)

      ! Defining zeros and weights for the Gauss-Hermite integration
      Allocate(XHERMI(1:number_Gauss),WHERMI(1:number_Gauss))
      Call HERMIT(number_Gauss,XHERMI,WHERMI)

      Allocate(PHERMI(1:number_HOwave+1),DHERMI(1:number_HOwave+1))

      Do i_Gauss=1,number_Gauss

         ! In the instruction below we have divided the Gauss weights by:
         !
         !             oscillator_scaling() = SQRT(MASS*HOMEG()/H_BAR**2)
         !
         ! This corresponds to expressing the one-dimensional integration
         ! measure in fermis. Integrals over the volume have therefore fermi**3
         ! as their dimension.
         WGTACT(i_Gauss,i_Cartesian)=WHERMI(i_Gauss)/scaling/oscillator_scaling(i_Cartesian)
         PNTACT(i_Gauss,i_Cartesian)=XHERMI(i_Gauss)/scaling

         XZER = PNTACT(i_Gauss,i_Cartesian)

         Call D_HERM(XZER,number_HOwave,PHERMI,DHERMI)

         ! The Hermite polynomials are multiplied by the dimensioned
         ! normalization constant:
         !
         !              (MASS*HOMEG()/H_BAR**2)**(1/4)
         !
         ! This corresponds to introducing the one-dimensional wave functions
         ! in units of fm**(-1/2). The corresponding densities in three
         ! dimensions are then in units of fm**(-3)

         Do i_wave=0,number_HOwave

            HERACT(i_wave,i_Gauss,i_Cartesian) = PHERMI(i_wave+1) &
                                               * Sqrt(oscillator_scaling(i_Cartesian)) &
                                               / HERFAC( i_wave )

            DHRACT(i_wave,i_Gauss,i_Cartesian) =(DHERMI(i_wave+1) - XZER*PHERMI(i_wave+1)) &
                                               * Sqrt(oscillator_scaling(i_Cartesian)) &
                                               / HERFAC( i_wave )

            DDHACT(i_wave,i_Gauss,i_Cartesian) =(XZER*XZER-2*i_wave-1) * PHERMI(i_wave+1) &
                                               * Sqrt(oscillator_scaling(i_Cartesian)) &
                                               / HERFAC( i_wave )

         End Do

         ! Derivatives of the "hermite * exponent" products (see comments in
         ! DEVHER) are multiplied by the appropriate powers of the dimensioned
         ! stretching factors "oscillator_scaling". In this way the first derivative has
         ! the correct dimension of fermi**(-3/2). Also the second derivative
         ! has the correct dimension of fermi**(-5/2).

         Do i_wave=0,number_HOwave
            DHRACT(i_wave,i_Gauss,i_Cartesian) = DHRACT(i_wave,i_Gauss,i_Cartesian) &
                                                              * oscillator_scaling(i_Cartesian)
            DDHACT(i_wave,i_Gauss,i_Cartesian) = DDHACT(i_wave,i_Gauss,i_Cartesian) &
                                                              * oscillator_scaling(i_Cartesian) &
                                                              * oscillator_scaling(i_Cartesian)
         End Do

      End Do

      Deallocate(XHERMI,WHERMI,PHERMI,DHERMI)

    End Subroutine build_HO_wavefunctions

    !---------------------------------------------------------------------!
    ! This subroutine calculates the Hermite polynomial normalization     !
    ! factors and stores them in common. it does the same as  norfac,     !
    ! which stors them in local array.                                    !
    !---------------------------------------------------------------------!
    Subroutine NORHER

      Integer(ipr) :: i

      If(Allocated(HERFAC)) Deallocate(HERFAC); Allocate(HERFAC(0:number_HOwave))

      HERFAC(0)=Sqrt(Sqrt(Pi))

      Do i=1,number_HOwave
         HERFAC(i)=HERFAC(i-1)*Sqrt(Real(2*i,Kind=pr))
      End Do

    End Subroutine NORHER

    !---------------------------------------------------------------------!
    !  Subroutine HERMIT calculates the zeroes X(i) of the NN-th order    !
    !  Hermite polynomial. The largest zero is stored in X(1). It also    !
    !  calculates the corresponding coefficients A(i) of the NN-th order  !
    !  Gauss-Hermite quadrature formula of degree 2*NN-1 .                !
    !  (Trivial modifications of the routine from Stroud book)            !
    !---------------------------------------------------------------------!
    Subroutine HERMIT(NN,X,A)

      Integer(ipr), INTENT(IN) :: NN
      Real(pr), Dimension(1:NN), INTENT(INOUT) :: X,A

      Integer(ipr) :: N1,N2,ierr,i,NI
      Real(pr) :: EEPS,CC,S,DDPN,PPN1,DPN,PN1,XT,XTT

      EEPS=precision_Hermite

      N1=NN-1; N2=(NN+1)/2

      CC=Sqrt(Pi)*FACTOR(NN-1)/(2.0_pr**N1)
      !S=(2.0_pr*Real(NN+1,Kind=pr))**0.16667_pr
      S=(2.0_pr*Real(NN+1,Kind=pr))**(1.0_pr/6.0_pr)

      Do i=1,N2

         If(i.Eq.1) XT=S**3 - 1.85575_pr/S
         If(i.Eq.2) XT=XT - 1.14_pr*NN**0.426_pr/XT
         If(i.Eq.3) XT=1.86_pr*XT - 0.86_pr*X(1)
         If(i.Eq.4) XT=1.91_pr*XT - 0.91_pr*X(2)
         If(i.Ge.5) XT=2.0_pr*XT - X(i-2)

         XTT=XT
         Call H_ROOT(XTT,NN,DDPN,PPN1,EEPS,ierr)

         DPN=DDPN; PN1=PPN1; XT=XTT

         X(i)=XT
         A(i)=CC/DPN/PN1

         NI=NN-i+1
         X(NI)=-XT
         A(NI)=A(i)

      End Do

    End Subroutine HERMIT

    !---------------------------------------------------------------------!
    ! This subroutine defines factorials                                  !
    !---------------------------------------------------------------------!
    Subroutine DEFFAC

      Integer(ipr) :: i

      FACTOR(0)=1.0_pr
      Do i=1,number_Gauss_max
         FACTOR(i)=FACTOR(i-1)*Real(i,Kind=pr)
      End Do

    End Subroutine DEFFAC

    !---------------------------------------------------------------------!
    ! Improves the approximated root X. In addition we also obtain        !
    !                                   PN1 = Value of H(N-1) at X        !
    ! (From the monograph by Stroud)                                      !
    !---------------------------------------------------------------------!

    Subroutine H_ROOT(X,NN,DPN,PN1,EPS,ierr)

      Integer(ipr), INTENT(IN) :: NN
      Integer(ipr), INTENT(INOUT) :: ierr
      Real(pr), INTENT(IN) :: EPS
      Real(pr), INTENT(INOUT) :: X,DPN,PN1

      Integer(ipr) :: ITERMX,ITER
      Real(pr) :: P,DP,D

      ITERMX=25

      Do ITER=1,ITERMX

         Call HRECUR(P,DP,PN1,X,NN)

         D=P/DP; X=X-D

         If(Abs(D).Le.EPS) Exit

      End Do

      DPN=DP

      If(ITER.Ge.ITERMX) Then
         ierr = 1
      Else
         ierr = 0
      End If

    End Subroutine H_ROOT

    !---------------------------------------------------------------------!
    ! Auxiliary recurrence routine (from the monograph by Stroud)         !
    !---------------------------------------------------------------------!
    Subroutine HRECUR(PN,DPN,PN1,X,NN)

      Integer(ipr), INTENT(IN) :: NN
      Real(pr), INTENT(INOUT) :: PN,DPN,PN1,X

      Integer(ipr) :: j
      Real(pr) :: P1,P,DP1,DP,FJ,FJ2,Q,DQ

      P1=1.0_pr; P=X
      DP1=0.0_pr; DP=1.0_pr

      Do j=2,NN

         FJ = Real(j,Kind=pr)
         FJ2 = 0.5_pr*(FJ - 1.0_pr)

         Q = X*P - FJ2*P1
         DQ = X*DP + P - FJ2*DP1

         P1=P; P=Q
         DP1=DP; DP=DQ

      End Do

      PN=P; DPN=DP; PN1=P1

    End Subroutine HRECUR

    !---------------------------------------------------------------------!
    ! Subroutine D_HERM calculates Hermite polynomials at point x up to   !
    ! n-th degree. It also calculates derivatives DHER(i) of these        !
    ! polynomials:                                                        !
    !                    DHER(I)=2*I*HER(I-1)                             !
    !---------------------------------------------------------------------!
    Subroutine D_HERM(X,N,HER,DHER)

      Integer(ipr), INTENT(IN) :: N
      Real(pr), INTENT(IN) :: X
      Real(pr), Dimension(1:N+1), INTENT(INOUT) :: HER,DHER

      Integer(ipr) :: i
      Real(pr) :: F,DF

      HER (1)=1.0_pr
      DHER(1)=0.0_pr

      If(N.Le.1) Return

      HER (2)=X + X
      DHER(2)=2.0_pr

      If(N.Le.2) Return

      Do i=2,N

         F=X*HER(i)-Real(i-1,Kind=pr)*HER(i-1)
         HER(i+1)=F+F

         DF=Real(i,Kind=pr)*HER(i)
         DHER(i+1)=DF+DF

      End Do

    End Subroutine D_HERM

    !---------------------------------------------------------------------!
    ! This subroutine calculates  the exponential factors:                !
    !                   EXP(-(ETA/SQRT(2))**2)                            !
    ! (eta denotes the gauss zeros).                                      !
    !---------------------------------------------------------------------!
    Subroutine build_HO_exponentials(i_Cartesian)

      Integer(ipr), INTENT(IN) :: i_Cartesian

      Integer(ipr) :: i

      Do i=1,number_Gauss
         EXPAUX(i,i_Cartesian)=Exp(-(PNTACT(i,i_Cartesian)**2))
      End Do

    End Subroutine build_HO_exponentials

    !---------------------------------------------------------------------!
    ! This subroutine copies the harmonic oscillator wave functions and   !
    ! their first and second derivatives to complex arrays.               !
    !---------------------------------------------------------------------!
    Subroutine COPHER(KARTEZ,NGAUSS)

      Integer(ipr), INTENT(IN) :: KARTEZ,NGAUSS

      Integer(ipr) :: IZEROS,NOSCIL,NMAIN2

      If(.Not.Allocated(CERMTS)) Allocate(CERMTS(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(CHRMTS)) Allocate(CHRMTS(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))
      If(.Not.Allocated(CDHRMT)) Allocate(CDHRMT(0:number_HOwave,1:number_Gauss_max,1:number_Cartesian))

      NMAIN2=2*NDOSCI+2

      Do IZEROS = 1,NGAUSS
         Do NOSCIL = 0,NMAIN2
            CERMTS(NOSCIL,IZEROS,KARTEZ) = HERACT(NOSCIL,IZEROS,KARTEZ)
            CHRMTS(NOSCIL,IZEROS,KARTEZ) = DHRACT(NOSCIL,IZEROS,KARTEZ)
            CDHRMT(NOSCIL,IZEROS,KARTEZ) = DDHACT(NOSCIL,IZEROS,KARTEZ)
         End Do
      End Do

     End Subroutine COPHER

    !---------------------------------------------------------------------!
    ! This subroutine calculates  the exponential factors:                !
    !                   EXP(-(ETA/SQRT(2))**2)                            !
    ! (eta denotes the gauss zeros).                                      !
    !---------------------------------------------------------------------!
    Subroutine SETBAS(NOSCIL,NLIMIT,ELIMIT)

      Integer(ipr), INTENT(IN) :: NOSCIL,NLIMIT
      Real(pr), INTENT(INOUT) :: ELIMIT

      Integer(ipr) :: MCOUNT,NX,NY,NZ,NUMCUT,MAXIMX,MAXIMY,MAXIMZ,MAXIXY,MAXIYZ,MAXIZX
      Integer(ipr), Dimension(NDBTOT) :: NBASIS

      Real(pr) :: EBASEX
      Real(pr), Dimension(NDBTOT) :: EBASIS

      !
      HOMSCA(1)=Sqrt(XMASSN*HBAROX)/H_BARC
      HOMSCA(2)=Sqrt(XMASSN*HBAROY)/H_BARC
      HOMSCA(3)=Sqrt(XMASSN*HBAROZ)/H_BARC

      ! Calculating the number of states in each of the 3 oscillators
      MCOUNT=0
      DO NX=0,NOSCIL
         DO NY=0,NOSCIL
            DO NZ=0,NOSCIL
               MCOUNT=MCOUNT+1
               NBASIS(MCOUNT)=MCOUNT
               EBASIS(MCOUNT)=HBAROX*(NX+0.5_pr)+HBAROY*(NY+0.5_pr)+HBAROZ*(NZ+0.5_pr)
            End Do
         End Do
      End Do

      NUMCUT=MCOUNT

      CALL sort_vector(EBASIS,NBASIS,MCOUNT)

      If (NLIMIT.LE.MCOUNT) Then
          ELIMIT=EBASIS(NLIMIT)+1.0D-5
      ELSE
          ELIMIT=EBASIS(MCOUNT)+1.0D-5
      End If

      INDICE(:,:,:)=0

      NUMCUT=0

      NXMAXX=0
      NYMAXX=0
      NZMAXX=0

      DO NX=0,NOSCIL
         DO NY=0,NOSCIL
            DO NZ=0,NOSCIL

               EBASEX=HBAROX*(NX+0.5_pr)+HBAROY*(NY+0.5_pr)+HBAROZ*(NZ+0.5_pr)

               If (EBASEX.LT.ELIMIT) Then

                   NUMCUT=NUMCUT+1

                   If (NUMCUT.LE.NDBASE) Then
                       NXVECT(NUMCUT)=NX
                       NYVECT(NUMCUT)=NY
                       NZVECT(NUMCUT)=NZ
                   End If

                   NXMAXX=Max(NXMAXX,NX)
                   NYMAXX=Max(NYMAXX,NY)
                   NZMAXX=Max(NZMAXX,NZ)

                   INDICE(NX,NY,NZ)=NUMCUT

               End If

            End Do
         End Do
      End Do

      LDBASE=NUMCUT

      ! X-DIRECTION
      DO NX=0,NXMAXX
         MAXIMY=0
         MAXIMZ=0
         DO NY=0,NYMAXX
            MAXIXY=0
            DO NZ=0,NZMAXX
               If (INDICE(NX,NY,NZ).NE.0) Then
                   MAXIMY=Max(MAXIMY,NY)
                   MAXIMZ=Max(MAXIMZ,NZ)
                   MAXIXY=Max(MAXIXY,NZ)
               End If
            End Do
            LAZOXY(NX,NY)=MAXIXY
         End Do
         LAYOFX(NX)=MAXIMY
         LAZOFX(NX)=MAXIMZ
      End Do

      ! Y-direction
      DO NY=0,NYMAXX
         MAXIMZ=0
         MAXIMX=0
         DO NZ=0,NZMAXX
            MAXIYZ=0
            DO NX=0,NXMAXX
               If (INDICE(NX,NY,NZ).NE.0) Then
                   MAXIMZ=Max(MAXIMZ,NZ)
                   MAXIMX=Max(MAXIMX,NX)
                   MAXIYZ=Max(MAXIYZ,NX)
               End If
            End Do
            LAXOYZ(NY,NZ)=MAXIYZ
         End Do
         LAZOFY(NY)=MAXIMZ
         LAXOFY(NY)=MAXIMX
      End Do

      !  Z-direction
      DO NZ=0,NZMAXX
         MAXIMX=0
         MAXIMY=0
         DO NX=0,NXMAXX
            MAXIZX=0
            DO NY=0,NYMAXX
               If (INDICE(NX,NY,NZ).NE.0) Then
                   MAXIMX=Max(MAXIMX,NX)
                   MAXIMY=Max(MAXIMY,NY)
                   MAXIZX=Max(MAXIZX,NY)
               End If
            End Do
            LAYOZX(NZ,NX)=MAXIZX
         End Do
         LAXOFZ(NZ)=MAXIMX
         LAYOFZ(NZ)=MAXIMY
      End Do

   End subroutine SETBAS

   !---------------------------------------------------------------------!
   !  This subroutine sorts a vector EL(1:N) by ascending order and      !
   !  returns the vector of indices in KL(1:N)                           !
   !---------------------------------------------------------------------!
   Subroutine sort_vector(EL,KL,N)

     Integer(ipr), INTENT(IN) :: N
     Integer(ipr), Dimension(1:N), INTENT(INOUT) :: KL
     Real(pr), Dimension(1:N), INTENT(INOUT) :: EL

     Integer(ipr) :: i,j,IEQUAL,IFILLL
     Integer(ipr), Allocatable :: KM(:),KWHER(:),KFROM(:)
     Real(pr), Allocatable :: EM(:)

     Allocate(KM(1:N),KWHER(1:N),KFROM(1:N))
     Allocate(EM(1:N))

     ! Rewriting input arrays to auxiliary ones
     KWHER(1:N)=1
     EM(1:N)=EL(1:N)
     KM(1:N)=KL(1:N)

     ! How many elements are smaller than the i-th one? In other words
     ! where should it go with respect to the offset position 1?
     Do J=1,N
        Do I=1,N
           If(EM(j).Lt.EM(i)) KWHER(i)=KWHER(i)+1
        End Do
     End Do

     ! Now we only have to resolve the case of identical elements
     ! From where the i-th element should come?
     KFROM(1:N)=0
     Do i=1,N
        KFROM(KWHER(i))=i
     End Do

     ! If some elements are equal there will be unfilled places
     ! in the KFROM array: Fill them up now.

   1 Continue

     Do i=1,N
        If(KFROM(i).Eq.0) Then
           IEQUAL=i-1
           IFILLL=i-1
           GO TO 2
        End IF
     End Do

     GO TO 3

   2 Do i=1,N

        If(KWHER(i).Eq.IEQUAL) Then
           KFROM(IFILLL)=i
           IFILLL=IFILLL+1
        End If

     End Do

     GO TO 1

   3 Continue

     ! Finally we may rewrite the input arrays in the right order
     Do i=1,N
        EL(i)=EM(KFROM(i))
        KL(i)=KM(KFROM(i))
     End Do

   End Subroutine sort_vector


 End Module HObasis