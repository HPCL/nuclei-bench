!-----------------------------------------------------------------------
!
!                MODULE DEFINING MPI COMMUNICATOR STRUCTURE
!
!-----------------------------------------------------------------------
#if(USE_MPI==1)
Module MPI_structure

   Use mpi
   Use hfodd_types

   Implicit None

   Integer, Public, Save :: worldGroup, groupMasters, groupSlaves       !< Groups
   Integer, Public, Save :: mastersCOMM, slavesCOMM, tribeCOMM          !< Communicators
   Integer, Public, Save :: numberMasters, numberHFODDproc, color       !< Number of tasks
   Integer, Public, Save :: tribeRank, slaveRank, masterRank, worldRank !< Ranks
   Integer, Public, Save :: tribeSize, slaveSize, masterSize, worldSize !< Sizes

Contains
   !=======================================================================
   !> This subroutine defines MPI environment for HFODD
   !=======================================================================
   Subroutine create_MPI_layout(cput_start)
      Real(pr), Intent(Inout) :: cput_start

      Integer(ipr) :: n_proc, n_tribe, n_color, intkey, mpi_err, ii, printTopology
      Integer(ipr), Dimension(:), Allocatable :: rankMasters
      Integer(ipr), Dimension(:,:), Allocatable :: rankTribe
      Integer(ipr), Dimension(:), Allocatable :: rankTribe1D
      Integer(ipr), Dimension(:), Allocatable :: all_tribeRank, all_color, all_worldRank

      Integer(ipr) :: NBLOCK,MBLOCK,NPGRID,MPGRID,KZHPEV
      COMMON &
             /SCALAPACK_INP/ NBLOCK,MBLOCK,NPGRID,MPGRID,KZHPEV

      ! Initializing the MPI environment and setting up MPI group structure
      Call mpi_init(mpi_err)
      Call mpi_comm_size(MPI_COMM_WORLD, worldSize, mpi_err)
      Call mpi_comm_rank(MPI_COMM_WORLD, worldRank, mpi_err)

      ! Record CPU starting times
      cput_start = MPI_Wtime()

      printTopology = 1
      numberHFODDproc = 1
      ! Attention: M_GRID and N_GRID below are passed through pre-processor
      !            not by standard NAMELI routine. This is because we need
       !            them even before we call NAMELI (or its likes)
#if(USE_MANYCORES==1)
      numberHFODDproc = M_GRID * N_GRID
      NPGRID = N_GRID
      MPGRID = M_GRID
#endif
      numberMasters = worldSize / numberHFODDproc

      ! Get handle on world group
      Call mpi_comm_group(MPI_COMM_WORLD, worldGroup, mpi_err)
      ! Define the rank of the masters
      Allocate(rankMasters(0:(numberMasters-1)))
      Do ii = 0, numberMasters - 1
         rankMasters(ii) = ii
      End Do
      ! Define a group of controllers
      Call mpi_group_incl(worldGroup, numberMasters, rankMasters, groupMasters, mpi_err)
      Call mpi_comm_create(MPI_COMM_WORLD, groupMasters, mastersCOMM, mpi_err)
      If(worldRank < numberMasters) Then
         Call mpi_comm_rank(mastersCOMM, masterRank, mpi_err)
         Call mpi_comm_size(mastersCOMM, masterSize, mpi_err)
      End If
      ! Everybody else is a worker
      If(numberHFODDproc > 1) Then
         Call mpi_group_excl(worldGroup, numberMasters, rankMasters, groupSlaves, mpi_err)
         Call mpi_comm_create(MPI_COMM_WORLD, groupSlaves, slavesCOMM, mpi_err)
         If(worldRank >= numberMasters) Then
            Call mpi_comm_rank(slavesCOMM, slaveRank, mpi_err)
            Call mpi_comm_size(slavesCOMM, slaveSize, mpi_err)
         Else
            slaveRank = 0
            slaveSize = 0
         End If
      End If

      ! Second grouping: define groups of workers associated with every controller. Since they are many such
      ! groups (as many as controllers), and they all have the same size, it is simpler to just split the
      ! communicator into many sub-communicators
      color = Mod(worldRank, numberMasters)
      intkey = worldRank
      Call mpi_comm_split(MPI_COMM_WORLD, color, intkey, tribeCOMM, mpi_err)
      Call mpi_comm_rank(tribeCOMM, tribeRank, mpi_err)
      Call mpi_comm_size(tribeCOMM, tribeSize, mpi_err)
      If(.Not.Allocated(all_tribeRank)) Allocate(all_tribeRank(0:(worldSize-1)),all_color(0:(worldSize-1)),all_worldRank(0:(worldSize-1)))
      Call MPI_Allgather(tribeRank, 1, MPI_INTEGER, all_tribeRank, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      Call MPI_Allgather(color, 1, MPI_INTEGER, all_color, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      Call MPI_Allgather(worldRank, 1, MPI_INTEGER, all_worldRank, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      If(numberMasters > 1) Then
         Allocate(rankTribe(0:(numberHFODDproc-1),0:(numberMasters-1)))
         Do n_proc = 0, worldSize - 1
            n_tribe = all_tribeRank(n_proc)
            n_color = all_color(n_proc)
            rankTribe(n_tribe, n_color) = all_worldRank(n_proc)
         End Do
         If(worldRank == 0) Then
            If(printTopology == 1) Then
               Do n_color = 0, numberMasters-1
                  Do n_tribe = 0, numberHFODDproc-1
                     Write(*,'("tribeRank = ",i4," n_color = ",i4, " worldRank = ",i4)') n_tribe,n_color,rankTribe(n_tribe, n_color)
                  End Do
               End Do
            End If
         End If
      Else
         Allocate(rankTribe1D(0:(numberHFODDproc-1)))
         Do n_proc = 0, worldSize - 1
            n_tribe = all_tribeRank(n_proc)
            rankTribe1D(n_tribe) = all_worldRank(n_proc)
         End Do
         If(worldRank == 0) Then
            If(printTopology == 1) Then
               Do n_tribe = 0, numberHFODDproc-1
                  Write(*,'("tribeRank = ",i4," worldRank = ",i4)') n_tribe,rankTribe1D(n_tribe)
               End Do
            End If
         End If
      End If
      Deallocate(all_tribeRank,all_color,all_worldRank)

   End Subroutine create_MPI_layout
   !=======================================================================
   !
   !=======================================================================
End Module MPI_structure
#endif

!-----------------------------------------------------------------------
!
!                MODULE DEFINING DENSHF_base: THE ROUTINE TO OPTIMIZE
!
!-----------------------------------------------------------------------
Module densities

   Use hfodd_types
   Use hfodd_sizes

   Implicit None

   Public DENSHF_base,allocate_data,deallocate_data

   Private

   Logical, Public, Save :: DOSYME=.True.

   Integer(ipr), Public, Save :: IHABOX=2,IHABOY=2,IHABOZ=2,IFTEMP=0

   Integer(ipr), Allocatable, Public, Save :: LDTOTS(:),LDSTAT(:),LDUPPE(:),LDTIMU(:)

   Complex(pr), Allocatable, Public, Save :: DE_RHO(:,:,:,:),DE_TAU(:,:,:,:),DE_LPR(:,:,:,:),DE_DIV(:,:,:,:),PD_RHO(:,:,:,:),PP_RHO(:,:,:,:)
   Complex(pr), Allocatable, Public, Save :: DE_SPI(:,:,:,:,:),DE_KIS(:,:,:,:,:),DE_GRR(:,:,:,:,:), &
                                             DE_LPS(:,:,:,:,:),DE_ROS(:,:,:,:,:),DE_ROC(:,:,:,:,:),DE_CUR(:,:,:,:,:)
   Complex(pr), Allocatable, Public, Save :: DE_SCU(:,:,:,:,:,:),DE_DES(:,:,:,:,:,:)
   Complex(pr), Allocatable, Public, Save :: WARIGH(:,:,:),WALEFT(:,:,:)
   Complex(pr), Allocatable, Public, Save :: SQUWAV(:,:,:,:)
   Complex(pr), Allocatable, Public, Save :: DCONTI(:,:,:)

Contains
   !=======================================================================
   !> Allocating all arrays. I simulate the use 'static' sizes, i.e., allocations are done with sizes that are
   !> greater than the ones actually used in the calculation. For instance, NDXHRM >= NXHERM, etc.
   !=======================================================================
   Subroutine allocate_data()
      Use HObasis, Only : NXVECT,NYVECT,NZVECT,INDICE,LAXOFY,LAXOFZ,LAYOFZ,LAYOFX,LAZOFX,LAZOFY,LAXOYZ,LAYOZX,LAZOXY,HOMSCA

      Allocate(NXVECT(NDBASE),NYVECT(NDBASE),NZVECT(NDBASE))
      Allocate(INDICE(0:NDXMAX,0:NDYMAX,0:NDZMAX))
      Allocate(LAXOFY(0:NDYMAX),LAXOFZ(0:NDZMAX))
      Allocate(LAYOFZ(0:NDZMAX),LAYOFX(0:NDXMAX))
      Allocate(LAZOFX(0:NDXMAX),LAZOFY(0:NDYMAX))
      Allocate(LAXOYZ(0:NDYMAX,0:NDZMAX))
      Allocate(LAYOZX(0:NDZMAX,0:NDXMAX))
      Allocate(LAZOXY(0:NDXMAX,0:NDYMAX))
      Allocate(HOMSCA(3))

      Allocate(DE_RHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(DE_TAU(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(DE_LPR(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(DE_DIV(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(PD_RHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(PP_RHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS))
      Allocate(DE_SPI(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_KIS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_GRR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_LPS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_ROS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_ROC(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_CUR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS))
      Allocate(DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1:NDKART,0:NDISOS))
      Allocate(DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1:NDKART,0:NDISOS))

      Allocate(WARIGH(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
      Allocate(WALEFT(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
      Allocate(LDTOTS(0:NDISOS),LDSTAT(0:NDISOS),LDUPPE(0:NDISOS),LDTIMU(0:NDISOS))

   End Subroutine allocate_data
   !=======================================================================
   !>  Deallocating all arrays
   !=======================================================================
   Subroutine deallocate_data()
      Use HObasis, Only : NXVECT,NYVECT,NZVECT,INDICE,LAXOFY,LAXOFZ,LAYOFZ,LAYOFX,LAZOFX,LAZOFY,LAXOYZ,LAYOZX,LAZOXY,HOMSCA

      Deallocate(NXVECT,NYVECT,NZVECT)
      Deallocate(INDICE)
      Deallocate(LAXOFY,LAXOFZ)
      Deallocate(LAYOFZ,LAYOFX)
      Deallocate(LAZOFX,LAZOFY)
      Deallocate(LAXOYZ,LAYOZX,LAZOXY)
      Deallocate(HOMSCA)
      Deallocate(DE_RHO,DE_TAU,DE_LPR,DE_DIV,PD_RHO,PP_RHO)
      Deallocate(DE_SPI,DE_KIS,DE_GRR,DE_LPS,DE_ROS,DE_ROC,DE_CUR)
      Deallocate(DE_SCU,DE_DES)
      Deallocate(WARIGH,WALEFT)
      Deallocate(LDTOTS,LDSTAT,LDUPPE,LDTIMU)

   End Subroutine deallocate_data
   !=======================================================================
   !>  Main routine
   !=======================================================================
   Subroutine DENSHF_base(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ, &
                          ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,ICHARG, &
                          IPNMIX,ITPNMX, &
                          NXMAXX,NAMEPN,PRINIT,IDEVAR,ISYMDE,IKERNE)
#if(USE_OPENMP==1)
      Use omp_lib
#endif
#if(USE_MPI==1 && USE_MANYCORES==1)
      Use MPI_structure
#endif
      Use HObasis, Only : LDBASE,CERMTS,CHRMTS,CDHRMT,INDICE,LAYOFX,LAZOXY

      Logical, Intent(In) :: PRINIT
      Integer(ipr), Intent(In) :: NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,IPAHFB, &
                                  MREVER,ICHARG,IPNMIX,ITPNMX,NXMAXX,IDEVAR,ISYMDE,IKERNE
      Character(Len=8), Intent(In) :: NAMEPN

      Logical :: DOSYME
      Integer(ipr) :: LREVER,IALLOC,I_FAIL,LTOTST,LSTATE,LUPPER,LTIMUP,LOFSET,LOFTRV,LUPTRV,LAYOFX_1,LAZOXY_1
      Integer(ipr) :: KX,KY,KZ,MU,NU,IX,IY,IZ,K,L,NX,NY,NZ,ISTATE,JSTATE,IBASE,ISPIN
      Integer(ipr), Dimension(0:NDKART) :: MUSYME
      Real(pr) :: U1OCCU,SIMFAC
      Real(pr), Dimension(1:2*NDSTAT) :: V1OCCU
      Complex(pr) :: ONEWAV,RMTS
      Complex(pr), Allocatable :: DENALL(:,:,:,:)
      Complex(pr), Dimension(:,:,:,:), Allocatable :: HRIGHE,HRIGDH,HRIGDD,HLEFHE,HLEFDH,HLEFDD
      Complex(pr), Dimension(:,:,:), Allocatable :: ARHEHE,ARHEDH,ARDHHE,ARHEDD,ARDDHE
      Complex(pr), Dimension(:,:,:), Allocatable :: ALHEHE,ALHEDH,ALDHHE,ALHEDD,ALDDHE
      Complex(pr), Dimension(-1:NDKART,-1:NDKART) :: DENAUX
      Complex(pr), Dimension(-1:NDKART,-1:NDKART) :: DEUPUP,DEDWDW,DEUPDW,DEDWUP
      Complex(pr), Dimension(-1:NDKART,-1:NDKART) :: PAUPUP,PADWDW,PAUPDW,PADWUP,PCUPUP,PCDWDW,PCUPDW,PCDWUP
      Complex(pr) :: DLUPUP,DLDWDW,DLUPDW,DLDWUP
      Complex(pr), Dimension(-1:NDKART,1:4*NDSTAT,0:NDSPIN) :: P,Q
      Complex(pr), Dimension(1:4*NDSTAT,0:NDSPIN) :: PLA,QLA
      Complex(pr) :: ZDOTU
#if(USE_MPI==1 && USE_MANYCORES==1)
      Integer(ipr) :: NSCALA,NVECTO,NTENSO
      Integer(ipr) :: mpi_err, buffer_size, i, j, block_size, block_mod, buffer, icount, KZ_all
      Integer(ipr), Allocatable :: block_vec(:), KZ_min(:), recvcnts(:), rdispls(:)
      Complex(pr), Allocatable :: bufsent(:), bufrecv(:)
#endif
      Logical :: DOSYMM = .False.

      !Call CPUTIM('DENSHF',1)

      LREVER=1

      IALLOC=0
      If(IDEVAR == 1 .And. IPAHFB /= 1) Then
         If(.NOT. Allocated(SQUWAV)) Then
            Allocate(SQUWAV(1:NDSTAT,NDXHRM,NDYHRM,NDZHRM),STAT=IALLOC)
            If(IALLOC /= 0) Call NOALLO('SQUWAV','DENSHF')
         End If
      End If
      If(.NOT. Allocated(DCONTI)) Then
         Allocate(DCONTI(1:NXHERM,1:NYHERM,1:NZHERM),STAT=IALLOC)
         If(IALLOC /= 0) Call NOALLO('DCONTI','DENSHF')
      End If
      Do KZ=1,NZHERM
         Do KY=1,NYHERM
            Do KX=1,NXHERM
               DCONTI(KX,KY,KZ)=C_ZERO
            End Do
         End Do
      End Do

      Do MU=0,NDKART
         If(IKERNE == 1) Then
            MUSYME(MU)=NDKART
         Else
            MUSYME(MU)=MU
         End If
      End Do

      ! Here we arrive with the wave functions that are already present in
      ! the WARIGH and WALEFT matrices
      LTOTST=LDTOTS(ICHARG)
      LSTATE=LDSTAT(ICHARG)
      LUPPER=LDUPPE(ICHARG)
      LTIMUP=LDTIMU(ICHARG)

      ! Symmetrization of densities is performed only for the most standard Skyrme calculations
      DOSYME=IPNMIX /= 1 .And. IKERNE /= 1 .And. .NOT. (ISYMDE == 1 .And. PRINIT) .And. DOSYMM

      ! Only half of the box in z-direction taken for ISIMTZ=1 or ISIGNY=1 or ISIQTY=1
      IHABOZ=2
      If((ISIMTZ /= 1 .And. ISIGNY /= 1 .And. ISIQTY /= 1) .Or. ( .NOT. DOSYME)) IHABOZ=1

      ! Only half of the box in y-direction taken for ISIMPY=1 or JSIMTY=1
      IHABOY=2
      If((ISIMPY /= 1 .And. JSIMTY /= 1) .Or. ( .NOT. DOSYME)) IHABOY=1

      !  Only half of the box in x-direction taken for ISIMTX=1
      IHABOX=2
      If(ISIMTX /= 1  .Or. ( .NOT. DOSYME)) IHABOX=1

      ! Loop over gauss-hermite points of z-integration
      PLA(:,:)=C_ZERO
      QLA(:,:)=C_ZERO
      P(:,:,:)=C_ZERO
      Q(:,:,:)=C_ZERO

      Allocate(ARHEHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),ARHEDH(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),&
               ARDHHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),ARHEDD(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),&
               ARDDHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN))
      Allocate(ALHEHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),ALHEDH(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),&
               ALDHHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),ALHEDD(0:NDXMAX,1:4*NDSTAT,0:NDSPIN),&
               ALDDHE(0:NDXMAX,1:4*NDSTAT,0:NDSPIN))
      If( .NOT. Allocated(HRIGHE)) Then
         I_FAIL=0
         Allocate(HRIGHE(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
         If(I_FAIL /= 0) Then
            Write(*,'("Error in Allocating HRIGHE - DENSHF")')
            Stop 'Error in Allocating HRIGHE'
         End If
         I_FAIL=0
         Allocate(HRIGDH(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
         If(I_FAIL /= 0) Then
            Write(*,'("Error in Allocating HRIGDH - DENSHF")')
            Stop 'Error in Allocating HRIGDH'
         End If
         I_FAIL=0
         Allocate(HRIGDD(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
         If(I_FAIL /= 0) Then
            Write(*,'("Error in Allocating HRIGDD - DENSHF")')
            Stop 'Error in Allocating HRIGDD'
         End If
         If(IKERNE == 1) Then
            I_FAIL=0
            Allocate(HLEFHE(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
            If(I_FAIL /= 0) Then
               Write(*,'("Error in Allocating HLEFHE - DENSHF")')
               Stop 'Error in Allocating HLEFHE'
            End If
            I_FAIL=0
            Allocate(HLEFDH(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
            If(I_FAIL /= 0) Then
               Write(*,'("Error in Allocating HLEFDH - DENSHF")')
               Stop 'Error in Allocating HLEFDH'
            End If
            I_FAIL=0
            Allocate(HLEFDD(0:NDXMAX,0:NDYMAX,1:4*NDSTAT,0:NDSPIN),STAT=I_FAIL)
            If(I_FAIL /= 0) Then
               Write(*,'("Error in Allocating HLEFDD - DENSHF")')
               Stop 'Error in Allocating HLEFDD'
            End If
         End If ! End If IKERNE
      End If ! End If ALLOCATED

      !=======================================================================
      !        LOOP OVER GAUSS-HERMITE POINTS OF Z-INTEGRATION
      !=======================================================================
#if(USE_MPI==1 && USE_MANYCORES==1)
      block_size = (NZHERM/IHABOZ)/numberHFODDproc
      block_mod = Mod(NZHERM/IHABOZ,numberHFODDproc)

      Allocate(block_vec(0:numberHFODDproc-1))
      If(block_mod > 0) Then
         buffer = block_mod
         Do i=0,numberHFODDproc-1
            If(buffer > 0) Then
               block_vec(i) = block_size + 1
               buffer=buffer-1
            Else
               block_vec(i) = block_size
            End If
         End Do
      Else
         Do i=0,numberHFODDproc-1
            block_vec(i) = block_size
         End Do
      End If

      Allocate(KZ_min(0:numberHFODDproc-1))
      KZ_min(0)=0
      Do i=1,numberHFODDproc-1
         KZ_min(i) = KZ_min(i-1) + block_vec(i-1)
      End Do

      Do IZ=KZ_min(tribeRank)+1,KZ_min(tribeRank)+block_vec(tribeRank)
#else
      Do IZ=1,NZHERM/IHABOZ
#endif
      ! The parallelization below seems to be the fastest. However, it it also the one that requires the
      ! largest amount of stack size (memory per  thread) because of the arrays HRIGHE, etc. For large
      ! bases and/or angular momentum calculations, it may be necessary to change the parallelization
      ! scheme so that multithreading takes place (i) first within the loop over NX (ii) Then also within
      ! the loop over IY
!$OMP PARALLEL DO &
!$OMP& DEFAULT(NONE) &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(NXMAXX,LTOTST,HRIGHE,HRIGDH,HRIGDD,HLEFHE,HLEFDH,HLEFDD,&
!$OMP&        LAYOFX,LAZOXY,CERMTS,CHRMTS,CDHRMT,IKERNE,INDICE,WARIGH,WALEFT,IZ) &
!$OMP& PRIVATE(NX,NY,NZ,ISPIN,IBASE,ISTATE,ONEWAV,RMTS,LAYOFX_1,LAZOXY_1)
         Do NX=0,NXMAXX
            LAYOFX_1=LAYOFX(NX)
            Do NY=0,LAYOFX_1
               Do ISTATE=1,LTOTST
                  Do ISPIN=0,NDSPIN

                     ONEWAV=C_ZERO
                     LAZOXY_1=LAZOXY(NX,NY)
                     Do NZ = 0,LAZOXY_1
                        IBASE = INDICE(NX,NY,NZ)
                        RMTS = CERMTS(NZ,IZ,3)
                        ONEWAV=ONEWAV+WARIGH(IBASE,ISTATE,ISPIN)*RMTS
                     End Do
                     HRIGHE(NX,NY,ISTATE,ISPIN)=ONEWAV
                     ONEWAV=C_ZERO
                     Do NZ = 0,LAZOXY_1
                        IBASE = INDICE(NX,NY,NZ)
                        RMTS = CHRMTS(NZ,IZ,3)
                        ONEWAV=ONEWAV+WARIGH(IBASE,ISTATE,ISPIN)*RMTS
                     End Do
                     HRIGDH(NX,NY,ISTATE,ISPIN)=ONEWAV
                     ONEWAV=C_ZERO
                     Do NZ = 0,LAZOXY_1
                        IBASE = INDICE(NX,NY,NZ)
                        RMTS = CDHRMT(NZ,IZ,3)
                        ONEWAV=ONEWAV+WARIGH(IBASE,ISTATE,ISPIN)*RMTS
                     End Do
                     HRIGDD(NX,NY,ISTATE,ISPIN)=ONEWAV

                     If(IKERNE == 1) Then
                        ONEWAV=C_ZERO
                        Do NZ = 0,LAZOXY_1
                           IBASE = INDICE(NX,NY,NZ)
                           RMTS = CERMTS(NZ,IZ,3)
                           ONEWAV=ONEWAV+WALEFT(IBASE,ISTATE,ISPIN)*RMTS
                        End Do
                        HLEFHE(NX,NY,ISTATE,ISPIN)=ONEWAV
                        ONEWAV=C_ZERO
                        Do NZ = 0,LAZOXY_1
                           IBASE = INDICE(NX,NY,NZ)
                           RMTS = CHRMTS(NZ,IZ,3)
                           ONEWAV=ONEWAV+WALEFT(IBASE,ISTATE,ISPIN)*RMTS
                        End Do
                        HLEFDH(NX,NY,ISTATE,ISPIN)=ONEWAV
                        ONEWAV=C_ZERO
                        Do NZ = 0,LAZOXY_1
                           IBASE = INDICE(NX,NY,NZ)
                           RMTS = CDHRMT(NZ,IZ,3)
                           ONEWAV=ONEWAV+WALEFT(IBASE,ISTATE,ISPIN)*RMTS
                        End Do
                        HLEFDD(NX,NY,ISTATE,ISPIN)=ONEWAV
                     End If

                  End Do !ISPIN
               End Do !ISTATE
            End Do !NY
         End Do !NX
!$OMP END PARALLEL DO
         !=======================================================================
         !        LOOP OVER GAUSS-HERMITE POINTS OF Y-INTEGRATION
         !=======================================================================
!$OMP PARALLEL DO &
!$OMP& DEFAULT(NONE) &
!$OMP& SCHEDULE(STATIC) &
!$OMP& SHARED(NXHERM,NYHERM,HRIGHE,HRIGDH,HRIGDD,HLEFHE,HLEFDH,HLEFDD, &
!$OMP&        IHABOX,IHABOY,NXMAXX,LSTATE,LTOTST,LAYOFX,LAZOXY, &
!$OMP&        CERMTS,CHRMTS,CDHRMT,IDEVAR,IPAHFB,IKERNE,V1OCCU, &
!$OMP&        SQUWAV,LUPPER,MUSYME,LREVER,ISIMPY, &
!$OMP&        LTIMUP,DE_RHO,DE_TAU,DE_LPR,DE_DIV,DE_SPI,DE_KIS,DE_LPS, &
!$OMP&        DE_GRR,DE_CUR,DE_ROS,DE_ROC,DE_SCU,DE_DES,PD_RHO,PP_RHO, &
!$OMP&        DCONTI,ITPNMX,MREVER,IZ,ICHARG) &
!$OMP& PRIVATE(IX,IY,NX,ISPIN,IBASE,ISTATE,MU,NU,Q,QLA,P,PLA,K,L, &
!$OMP&         ARHEHE,ARHEDH,ARDHHE,ARHEDD,ARDDHE,ALHEHE,ALHEDH,ALDHHE, &
!$OMP&         ALHEDD,ALDDHE,JSTATE,U1OCCU,DLUPUP,DLDWDW,DLUPDW,DLDWUP, &
!$OMP&         DEUPUP,DEDWDW,DEUPDW,DEDWUP,DENAUX,LOFSET,PAUPUP,PADWDW, &
!$OMP&         PAUPDW,PADWUP,PCUPUP,PCDWDW,PCUPDW,PCDWUP,LOFTRV,LUPTRV, &
!$OMP&         SIMFAC,LAYOFX_1)
         Do IY=1,NYHERM/IHABOY

            Do ISPIN=0,NDSPIN
               Do ISTATE=1,LTOTST
                  Do NX=0,NXMAXX
                     LAYOFX_1 = LAYOFX(NX) + 1
                     ARHEHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HRIGHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                     ARHEDH(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HRIGDH(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                     ARDHHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HRIGHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CHRMTS(0,IY,2),1)
                     ARHEDD(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HRIGDD(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                     ARDDHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HRIGHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CDHRMT(0,IY,2),1)
                     If(IKERNE == 1) Then
                        ALHEHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HLEFHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                        ALHEDH(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HLEFDH(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                        ALDHHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HLEFHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CHRMTS(0,IY,2),1)
                        ALHEDD(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HLEFDD(NX,0,ISTATE,ISPIN),NDXMAX+1,CERMTS(0,IY,2),1)
                        ALDDHE(NX,ISTATE,ISPIN)=ZDOTU(LAYOFX_1,HLEFHE(NX,0,ISTATE,ISPIN),NDXMAX+1,CDHRMT(0,IY,2),1)
                     End If
                  End Do ! end loop NX
               End Do ! end loop ISTATE
            End Do ! end loop ISPIN

            !=======================================================================
            !           LOOP OVER GAUSS-HERMITE POINTS OF X-INTEGRATION
            !=======================================================================

            Do IX=1,NXHERM/IHABOX

               Do ISPIN=0,NDSPIN

                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARHEHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,Q(0,1,ISPIN), NDKART+2)
                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARHEHE(0,1,ISPIN),NDXMAX+1,CHRMTS(0,IX,1),1,C_ZERO,Q(1,1,ISPIN), NDKART+2)
                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARDHHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,Q(2,1,ISPIN), NDKART+2)
                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARHEDH(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,Q(3,1,ISPIN), NDKART+2)

                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARHEHE(0,1,ISPIN),NDXMAX+1,CDHRMT(0,IX,1),1,C_ZERO,QLA(1,ISPIN), 1)
                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARDDHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_UNIT,QLA(1,ISPIN), 1)
                  Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                             ARHEDD(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_UNIT,QLA(1,ISPIN), 1)

                  If(IKERNE == 1) Then

                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALHEHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,P(0,1,ISPIN), NDKART+2)
                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALHEHE(0,1,ISPIN),NDXMAX+1,CHRMTS(0,IX,1),1,C_ZERO,P(1,1,ISPIN), NDKART+2)
                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALDHHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,P(2,1,ISPIN), NDKART+2)
                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT,&
                                ALHEDH(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_ZERO,P(3,1,ISPIN), NDKART+2)

                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALHEHE(0,1,ISPIN),NDXMAX+1,CDHRMT(0,IX,1),1,C_ZERO,PLA(1,ISPIN), 1)
                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALDDHE(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_UNIT,PLA(1,ISPIN), 1)
                     Call ZGEMV('T',NXMAXX+1,LTOTST,C_UNIT, &
                                ALHEDD(0,1,ISPIN),NDXMAX+1,CERMTS(0,IX,1),1,C_UNIT,PLA(1,ISPIN), 1)

                  End If ! End If IKERNE

               End Do ! end loop ISPIN

               ! At this point:
               !  - Q(0,ISTATE,ISPIN) contains the ISTATE-th single-particle "right" wave function, or upper/lower
               !                      quasi-particle wave function, in one of the spin blocks, as described above.
               !  - Q(K,ISTATE,ISPIN) (K=1,2,3) Contains its x,y,z derivative
               !  - QLA(ISTATE,ISPIN) contains its Laplacian
               !
               !  - P(0,ISTATE,ISPIN) contains the istate-th single-particle "left" wave function, or upper/lower
               !                      quasi-particle wave function, in one of the spin blocks, as described above.
               !  - P(K,ISTATE,ISPIN) (K=1,2,3) Contains its x,y,z derivative
               !  - PLA(ISTATE,ISPIN) Contains its laplacian
               ! They are calculated at the (IX,IY,IZ)-th point of the integration mesh

               ! For IKERNE=0 Q( . . . ) is stored in P( . . . )
               If(IKERNE /= 1) Then
                  Call ZCOPY((NDKART+1)*4*NDSTAT*(NDSPIN+1),Q(0,1,0),1,P(0,1,0),1)
                  Call ZCOPY(4*NDSTAT*(NDSPIN+1),QLA,1,PLA,1)
               End If

               ! Complex conjugate of Q( . . . ) is stored in Q
               Do ISPIN=0,NDSPIN
                  Do ISTATE=1,LTOTST
                     Do MU=0,NDKART
                        Q(MU,ISTATE,ISPIN)=Conjg(Q(MU,ISTATE,ISPIN))
                     End Do
                  End Do
               End Do

               Do ISPIN=0,NDSPIN
                  Do ISTATE=1,LTOTST
                     QLA(ISTATE,ISPIN)=Conjg( QLA(ISTATE,ISPIN))
                  End Do
               End Do

               ! Storing squares of wave functions for the calculation of pairing matrix elements
               If(IDEVAR == 1 .And. IPAHFB < 1) Then
                  If(LSTATE > NDSTAT) STOP ' LSTATE > NDSTAT IN DENSHF'
                  JSTATE=LSTATE
                  Do ISTATE=1,LSTATE
                     JSTATE=JSTATE+1
                     U1OCCU=Sqrt(1.0-V1OCCU(ISTATE)**2)
                     SQUWAV(ISTATE,IX,IY,IZ) = 2.0*Real(P(0,ISTATE,1)*Q(0,ISTATE,1) + P(0,ISTATE,0)*Q(0,ISTATE,0))
                     ! Here the wave functions are multiplied by the occupation probabilties
                     Do ISPIN=0,NDSPIN
                        Call ZDSCAL(NDKART+1,V1OCCU(ISTATE),P(0,ISTATE,ISPIN),1)
                        Call ZDSCAL(NDKART+1,V1OCCU(ISTATE),Q(0,ISTATE,ISPIN),1)
                        PLA(ISTATE,ISPIN) = PLA(ISTATE,ISPIN)*V1OCCU(ISTATE)
                        QLA(ISTATE,ISPIN) = QLA(ISTATE,ISPIN)*V1OCCU(ISTATE)
                     End Do ! end loop ISPIN
                  End Do ! end loop ISTATE
               End If ! End If IDEVAR

               !=======================================================================
               !              SUMMING UP DENSITIES IN THE P-H CHANNEL
               !=======================================================================

               ! Summing up contributions to the auxiliary p-h densities
               DLUPUP = C_ZERO
               DLDWDW = C_ZERO
               DLUPDW = C_ZERO
               DLDWUP = C_ZERO

               DLUPUP = (ZDOTU(LUPPER,PLA(1,1),1,Q(0,1,1),NDKART+2) + ZDOTU(LUPPER,QLA(1,1),1,P(0,1,1),NDKART+2))/2
               DLDWDW = (ZDOTU(LUPPER,PLA(1,0),1,Q(0,1,0),NDKART+2) + ZDOTU(LUPPER,QLA(1,0),1,P(0,1,0),NDKART+2))/2
               DLUPDW = (ZDOTU(LUPPER,PLA(1,1),1,Q(0,1,0),NDKART+2) + ZDOTU(LUPPER,QLA(1,0),1,P(0,1,1),NDKART+2))/2
               DLDWUP = (ZDOTU(LUPPER,PLA(1,0),1,Q(0,1,1),NDKART+2) + ZDOTU(LUPPER,QLA(1,1),1,P(0,1,0),NDKART+2))/2

               ! p-h densities

               ! Spin up-up and down-down p-h densities
               Call ZGEMM('N','T',NDKART+1,NDKART+1,LUPPER,C_UNIT, &
                          P(0,1,1),NDKART+2,Q(0,1,1),NDKART+2,C_ZERO,DEUPUP(0,0), NDKART+2)
               Call ZGEMM('N','T',NDKART+1,NDKART+1,LUPPER,C_UNIT, &
                          P(0,1,0),NDKART+2,Q(0,1,0),NDKART+2, C_ZERO,DEDWDW(0,0), NDKART+2)

               ! Spin up-down p-h densities
               Call ZGEMM('N','T',NDKART+1,NDKART+1,LUPPER,C_UNIT, &
                          P(0,1,1),NDKART+2,Q(0,1,0),NDKART+2,C_ZERO,DEUPDW(0,0),NDKART+2)

               ! Spin down-up p-h densities
               If(IKERNE == 1) Then
                  Call ZGEMM('N','T',NDKART+1,NDKART+1,LUPPER, &
                             C_UNIT,P(0,1,0),NDKART+2,Q(0,1,1),NDKART+2,C_ZERO,DEDWUP(0,0),NDKART+2)
               Else
                  Do MU=0,NDKART
                     Do NU=0,NDKART
                        DEDWUP(MU,NU)=Conjg(DEUPDW(NU,MU))
                     End Do
                  End Do
               End If

               ! Scalar p.h. densities
               DE_RHO(IX,IY,IZ,ITPNMX)=DE_RHO(IX,IY,IZ,ITPNMX)+DEUPUP(0,0)+DEDWDW(0,0)

               DE_TAU(IX,IY,IZ,ITPNMX)=DE_TAU(IX,IY,IZ,ITPNMX)+DEUPUP(1,1)+DEDWDW(1,1) &
                                                              +DEUPUP(2,2)+DEDWDW(2,2) &
                                                              +DEUPUP(3,3)+DEDWDW(3,3)

               DE_LPR(IX,IY,IZ,ITPNMX)=DE_LPR(IX,IY,IZ,ITPNMX)+2*(DLUPUP     +DLDWDW &
                                                                 +DEUPUP(1,1)+DEDWDW(1,1) &
                                                                 +DEUPUP(2,2)+DEDWDW(2,2) &
                                                                 +DEUPUP(3,3)+DEDWDW(3,3))

               DE_DIV(IX,IY,IZ,ITPNMX)=DE_DIV(IX,IY,IZ,ITPNMX)+((DEUPDW(2,3)-DEUPDW(3,2)) &
                                                               +(DEDWUP(2,3)-DEDWUP(3,2)))*UNIT_I &
                                                              +((DEUPDW(3,1)-DEUPDW(1,3)) &
                                                              - (DEDWUP(3,1)-DEDWUP(1,3)))*(-1) &
                                                              +((DEUPUP(1,2)-DEUPUP(2,1)) &
                                                               -(DEDWDW(1,2)-DEDWDW(2,1)))*UNIT_I

               ! Vector p.h. densities
               DE_SPI(IX,IY,IZ,1,ITPNMX)=DE_SPI(IX,IY,IZ,1,ITPNMX)+(DEUPDW(0,0)+DEDWUP(0,0))
               DE_SPI(IX,IY,IZ,2,ITPNMX)=DE_SPI(IX,IY,IZ,2,ITPNMX)+(DEUPDW(0,0)-DEDWUP(0,0))*UNIT_I
               DE_SPI(IX,IY,IZ,3,ITPNMX)=DE_SPI(IX,IY,IZ,3,ITPNMX)+(DEUPUP(0,0)-DEDWDW(0,0))

               DE_KIS(IX,IY,IZ,1,ITPNMX)=DE_KIS(IX,IY,IZ,1,ITPNMX)+(DEUPDW(1,1)+DEDWUP(1,1) &
                                                                   +DEUPDW(2,2)+DEDWUP(2,2) &
                                                                   +DEUPDW(3,3)+DEDWUP(3,3))
               DE_KIS(IX,IY,IZ,2,ITPNMX)=DE_KIS(IX,IY,IZ,2,ITPNMX)+(DEUPDW(1,1)-DEDWUP(1,1) &
                                                                   +DEUPDW(2,2)-DEDWUP(2,2) &
                                                                   +DEUPDW(3,3)-DEDWUP(3,3))*UNIT_I
               DE_KIS(IX,IY,IZ,3,ITPNMX)=DE_KIS(IX,IY,IZ,3,ITPNMX)+(DEUPUP(1,1)-DEDWDW(1,1) &
                                                                   +DEUPUP(2,2)-DEDWDW(2,2) &
                                                                   +DEUPUP(3,3)-DEDWDW(3,3))

               DE_LPS(IX,IY,IZ,1,ITPNMX)=DE_LPS(IX,IY,IZ,1,ITPNMX)+2*(DLUPDW     +DLDWUP &
                                                                     +DEUPDW(1,1)+DEDWUP(1,1) &
                                                                     +DEUPDW(2,2)+DEDWUP(2,2) &
                                                                     +DEUPDW(3,3)+DEDWUP(3,3))

               DE_LPS(IX,IY,IZ,2,ITPNMX)=DE_LPS(IX,IY,IZ,2,ITPNMX)+2*(DLUPDW     -DLDWUP &
                                                                     +DEUPDW(1,1)-DEDWUP(1,1) &
                                                                     +DEUPDW(2,2)-DEDWUP(2,2) &
                                                                     +DEUPDW(3,3)-DEDWUP(3,3))*UNIT_I

               DE_LPS(IX,IY,IZ,3,ITPNMX)=DE_LPS(IX,IY,IZ,3,ITPNMX)+2*(DLUPUP     -DLDWDW &
                                                                     +DEUPUP(1,1)-DEDWDW(1,1) &
                                                                     +DEUPUP(2,2)-DEDWDW(2,2) &
                                                                     +DEUPUP(3,3)-DEDWDW(3,3))

               Do K=1,NDKART

                  DE_GRR(IX,IY,IZ,K,ITPNMX)=DE_GRR(IX,IY,IZ,K,ITPNMX)+(DEUPUP(K,0)+DEUPUP(0,K) &
                                                                      +DEDWDW(K,0)+DEDWDW(0,K))

                  DE_CUR(IX,IY,IZ,K,ITPNMX)=DE_CUR(IX,IY,IZ,K,ITPNMX)+(DEUPUP(K,0)-DEUPUP(0,K) &
                                                                      +DEDWDW(K,0)-DEDWDW(0,K))/2/UNIT_I

                  DENAUX(K,1)=((DEUPDW(K,0)+DEUPDW(0,K))+(DEDWUP(K,0)+DEDWUP(0,K)))
                  DENAUX(K,2)=((DEUPDW(K,0)+DEUPDW(0,K))-(DEDWUP(K,0)+DEDWUP(0,K)))*UNIT_I
                  DENAUX(K,3)=((DEUPUP(K,0)+DEUPUP(0,K))-(DEDWDW(K,0)+DEDWDW(0,K)))

               End Do

               DE_ROS(IX,IY,IZ,1,ITPNMX)=DE_ROS(IX,IY,IZ,1,ITPNMX)+DENAUX(2,3)-DENAUX(3,2)
               DE_ROS(IX,IY,IZ,2,ITPNMX)=DE_ROS(IX,IY,IZ,2,ITPNMX)+DENAUX(3,1)-DENAUX(1,3)
               DE_ROS(IX,IY,IZ,3,ITPNMX)=DE_ROS(IX,IY,IZ,3,ITPNMX)+DENAUX(1,2)-DENAUX(2,1)

               DE_ROC(IX,IY,IZ,1,ITPNMX)=DE_ROC(IX,IY,IZ,1,ITPNMX)+(DEUPUP(2,3)-DEUPUP(3,2) &
                                                                   +DEDWDW(2,3)-DEDWDW(3,2))*UNIT_I

               DE_ROC(IX,IY,IZ,2,ITPNMX)=DE_ROC(IX,IY,IZ,2,ITPNMX)+(DEUPUP(3,1)-DEUPUP(1,3) &
                                                                   +DEDWDW(3,1)-DEDWDW(1,3))*UNIT_I

               DE_ROC(IX,IY,IZ,3,ITPNMX)=DE_ROC(IX,IY,IZ,3,ITPNMX)+(DEUPUP(1,2)-DEUPUP(2,1) &
                                                                   +DEDWDW(1,2)-DEDWDW(2,1))*UNIT_I

               ! Tensor p.h. densities
               Do K=1,NDKART
                  DE_SCU(IX,IY,IZ,K,1,ITPNMX)=DE_SCU(IX,IY,IZ,K,1,ITPNMX) &
                                           +((DEUPDW(K,0)-DEUPDW(0,K))+(DEDWUP(K,0)-DEDWUP(0,K)))/2/UNIT_I
                  DE_SCU(IX,IY,IZ,K,2,ITPNMX)=DE_SCU(IX,IY,IZ,K,2,ITPNMX) &
                                           +((DEUPDW(K,0)-DEUPDW(0,K))-(DEDWUP(K,0)-DEDWUP(0,K)))/2
                  DE_SCU(IX,IY,IZ,K,3,ITPNMX)=DE_SCU(IX,IY,IZ,K,3,ITPNMX) &
                                           +((DEUPUP(K,0)-DEUPUP(0,K))-(DEDWDW(K,0)-DEDWDW(0,K)))/2/UNIT_I
                  Do L=1,NDKART
                     DE_DES(IX,IY,IZ,K,L,ITPNMX)=DE_DES(IX,IY,IZ,K,L,ITPNMX)+DENAUX(K,L)
                  End Do
               End Do

               !=======================================================================
               !             SUMMING UP DENSITIES IN THE P-P CHANNEL
               !=======================================================================

               If(IPAHFB == 1) Then

                  LOFSET=0
                  LOFTRV=LTIMUP
                  LUPTRV=LSTATE
                  If(LREVER == 0) Then
                     LOFSET=LTIMUP
                     LOFTRV=0
                  End If

                  ! "Conjugate" pairing density, < L|a+a+|R > <--> kappa^{01}, PD_RHO
                  ! Time-up contribution
                  !  Additional explanation for PNP in the canonical basis
                  !   * Canonical wave functions for i > LUPPER contain (2\sigma') \phi^{*}_{n}(r,-\sigma')
                  !   * P(0,:,:) contains the canonical wave function i <= LUPPER, hence v_n \phi_n(r,\sigma)
                  !   * Q(0,:,:) contains the complex conjugate of the canonical wave function i > LUPPER, hence
                  !     u_n /(u_n^2 + z^2 v_n^2) (2\sigma') \phi_n(r,-\sigma')
                  !  Therefore: (u_n v_n) /(u_n^2 + z^2 v_n^2) (2\sigma') \phi_n(r,\sigma)\phi^_n(r,-\sigma')
                  Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                             P(0,LOFSET+1,1),NDKART+2,Q(0,LUPPER+1,1),NDKART+2,C_ZERO,PAUPUP(0,0),NDKART+2)
                  Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                             P(0,LOFSET+1,0),NDKART+2,Q(0,LUPPER+1,0),NDKART+2,C_ZERO,PADWDW(0,0),NDKART+2)
                  Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                             P(0,LOFSET+1,1),NDKART+2,Q(0,LUPPER+1,0),NDKART+2,C_ZERO,PAUPDW(0,0),NDKART+2)
                  Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                             P(0,LOFSET+1,0),NDKART+2,Q(0,LUPPER+1,1),NDKART+2,C_ZERO,PADWUP(0,0),NDKART+2)

                  ! Time-down contribution
                  If(IKERNE == 1 .Or. ISIMPY == 0) Then
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LOFTRV+1,1),NDKART+2,Q(0,LUPTRV+1,1),NDKART+2,C_UNIT,PAUPUP(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LOFTRV+1,0),NDKART+2,Q(0,LUPTRV+1,0),NDKART+2,C_UNIT,PADWDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LOFTRV+1,1),NDKART+2,Q(0,LUPTRV+1,0),NDKART+2,C_UNIT,PAUPDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LOFTRV+1,0),NDKART+2,Q(0,LUPTRV+1,1),NDKART+2,C_UNIT,PADWUP(0,0),NDKART+2)
                  End If

                  ! "Regular" pairing density, < L|aa|R > <--> kappa^{10}, PP_RHO
                  If(IKERNE == 1 .And. MREVER /= 1) STOP "IKERNE=1 WITH IPAHFB=1 REQUIRES MREVER=1 IN DENSHF"

                  !  Additional explanation for PNP in the canonical basis
                  !   * Canonical wave functions for i > LUPPER contain (2\sigma') \phi^{*}_{n}(r,-\sigma')
                  !   * P(0,:,:) contains the canonical wave function i > LUPPER, hence u_n (2\sigma) \phi^{*}_n(r,-\sigma)
                  !   * Q(0,:,:) contains the complex conjugate of the canonical wave function i <= LUPPER, hence
                  !     (z^2 v_n) /(u_n^2 + z^2 v_n^2) \phi^{*}_{n}(r,\sigma')
                  !  Therefore: (u_n v_n z^2) /(u_n^2 + z^2 v_n^2) (2\sigma) \phi^{*}_{n}(r,-\sigma)\phi^{*}_{n}(r,\sigma')
                  If(IKERNE == 1) Then
                     ! Time-up contribution
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                                P(0,LUPPER+1,1),NDKART+2,Q(0,LOFSET+1,1),NDKART+2,C_ZERO,PCUPUP(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                                P(0,LUPPER+1,0),NDKART+2,Q(0,LOFSET+1,0),NDKART+2,C_ZERO,PCDWDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                                P(0,LUPPER+1,1),NDKART+2,Q(0,LOFSET+1,0),NDKART+2,C_ZERO,PCUPDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LSTATE-LUPPER,C_UNIT, &
                                P(0,LUPPER+1,0),NDKART+2,Q(0,LOFSET+1,1),NDKART+2,C_ZERO,PCDWUP(0,0),NDKART+2)
                     ! Time-down contribution
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LUPTRV+1,1),NDKART+2,Q(0,LOFTRV+1,1),NDKART+2,C_UNIT,PCUPUP(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LUPTRV+1,0),NDKART+2,Q(0,LOFTRV+1,0),NDKART+2,C_UNIT,PCDWDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LUPTRV+1,1),NDKART+2,Q(0,LOFTRV+1,0),NDKART+2,C_UNIT,PCUPDW(0,0),NDKART+2)
                     Call ZGEMM('N','T',NDKART+1,NDKART+1,LTOTST-LSTATE,C_UNIT, &
                                P(0,LUPTRV+1,0),NDKART+2,Q(0,LOFTRV+1,1),NDKART+2,C_UNIT,PCDWUP(0,0),NDKART+2)
                  Else
                     PCUPUP(:,:) = Conjg(PAUPUP(:,:))
                     PCDWDW(:,:) = Conjg(PADWDW(:,:))
                     PCUPDW(:,:) = Conjg(PAUPDW(:,:))
                     PCDWUP(:,:) = Conjg(PADWUP(:,:))
                  End If

                  ! For the conserved simplex (ISIMPY=1) only half of the quasiparticles contribute to
                  ! the p-p densities, and therefore the densities below are then multiplied by 2.
                  ! For pairing functionals restricted to time-even terms, only the real part of the
                  ! pairing density is kept.

                  ! Multiply by two except if kernel calculation, in which cases we used the whole
                  ! basis to calculate the pair densities
                                                   SIMFAC=2.0
                  If(IKERNE == 1 .Or. ISIMPY == 0) SIMFAC=1.0

                  ! MREVER==1 for broken T or kernel
                  If(MREVER == 1) Then
                     PD_RHO(IX,IY,IZ,ITPNMX) = PD_RHO(IX,IY,IZ,ITPNMX) + SIMFAC*(PAUPUP(0,0)+PADWDW(0,0))
                     PP_RHO(IX,IY,IZ,ITPNMX) = PP_RHO(IX,IY,IZ,ITPNMX) + SIMFAC*(PCUPUP(0,0)+PCDWDW(0,0))
                  ! MREVER==0 otherwise. Pairing density is then real
                  Else
                     PD_RHO(IX,IY,IZ,ITPNMX) = PD_RHO(IX,IY,IZ,ITPNMX) + SIMFAC*REAL(PAUPUP(0,0)+PADWDW(0,0))
                     PP_RHO(IX,IY,IZ,ITPNMX) = PP_RHO(IX,IY,IZ,ITPNMX) + SIMFAC*REAL(PCUPUP(0,0)+PCDWDW(0,0))
                  End If

               End If

               !=======================================================================
               !              ENDING THE DO-LOOPS OVER THE CARTESIAN POINTS
               !=======================================================================

            End Do ! END OF LOOP OVER IX
         End Do ! END OF LOOP OVER IY
!$OMP END PARALLEL DO
      End Do ! END OF LOOP OVER IZ

      Deallocate(ARHEHE,ARHEDH,ARDHHE,ARHEDD,ARDDHE)
      Deallocate(ALHEHE,ALHEDH,ALDHHE,ALHEDD,ALDDHE)

      Deallocate(HRIGHE,HRIGDH,HRIGDD,STAT=I_FAIL)
      If(I_FAIL /= 0) Then
         Write(*,'("Error in DEALLOCATING arrays - DENSHF")')
         Stop 'Error in DEALLOCATION'
      End If
      If(IKERNE == 1) Then
         Deallocate(HLEFHE,HLEFDH,HLEFDD,STAT=I_FAIL)
         If(I_FAIL /= 0) Then
            Write(*,'("Error in DEALLOCATING arrays - DENSHF")')
            Stop 'Error in DEALLOCATION'
         End If
      End If

#if(USE_MPI==1 && USE_MANYCORES==1)

      NSCALA=6
      NVECTO=7
      NTENSO=2

      Allocate(recvcnts(0:numberHFODDproc-1))
      Do i=0,numberHFODDproc-1
         recvcnts(i) = block_vec(i) * (NYHERM/IHABOY) * (NXHERM/IHABOX) * (NSCALA + NVECTO*3 + NTENSO*9)
      End Do

      buffer_size =  recvcnts(tribeRank)
      Allocate(bufsent(buffer_size))

      icount=0
      ! scalar densities
      Do IZ=KZ_min(tribeRank)+1,KZ_min(tribeRank)+block_vec(tribeRank)
         Do IY=1,NYHERM/IHABOY
            Do IX=1,NXHERM/IHABOX
               icount=icount+1
               bufsent(icount)=DE_RHO(IX,IY,IZ,ITPNMX)
               icount=icount+1
               bufsent(icount)=DE_TAU(IX,IY,IZ,ITPNMX)
               icount=icount+1
               bufsent(icount)=DE_LPR(IX,IY,IZ,ITPNMX)
               icount=icount+1
               bufsent(icount)=DE_DIV(IX,IY,IZ,ITPNMX)
               icount=icount+1
               bufsent(icount)=PD_RHO(IX,IY,IZ,ITPNMX)
               icount=icount+1
               bufsent(icount)=PP_RHO(IX,IY,IZ,ITPNMX)
            End Do
         End Do
      End Do
      ! vector densities
      Do K=1,NDKART
         Do IZ=KZ_min(tribeRank)+1,KZ_min(tribeRank)+block_vec(tribeRank)
            Do IY=1,NYHERM/IHABOY
               Do IX=1,NXHERM/IHABOX
                  icount=icount+1
                  bufsent(icount)=DE_SPI(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_KIS(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_LPS(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_ROS(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_ROC(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_GRR(IX,IY,IZ,K,ITPNMX)
                  icount=icount+1
                  bufsent(icount)=DE_CUR(IX,IY,IZ,K,ITPNMX)
               End Do
            End Do
         End Do
      End Do
      ! tensor densities
      Do L=1,NDKART
         Do K=1,NDKART
            Do IZ=KZ_min(tribeRank)+1,KZ_min(tribeRank)+block_vec(tribeRank)
               Do IY=1,NYHERM/IHABOY
                  Do IX=1,NXHERM/IHABOX
                     icount=icount+1
                     bufsent(icount)=DE_SCU(IX,IY,IZ,K,L,ITPNMX)
                     icount=icount+1
                     bufsent(icount)=DE_DES(IX,IY,IZ,K,L,ITPNMX)
                  End Do
               End Do
            End Do
          End Do
      End Do

      Allocate(rdispls(0:numberHFODDproc-1))
      rdispls(0)=0
      Do i=1,numberHFODDproc-1
         rdispls(i) = rdispls(i-1) + recvcnts(i-1)
      End Do

      KZ_all = (NZHERM/IHABOZ) * (NYHERM/IHABOY) * (NXHERM/IHABOX) * (NSCALA + NVECTO*3 + NTENSO*9)
      If(KZ_all /= SUM(recvcnts)) Then
         Write(*,'("KZ_all=",i15," SUM(recvcnts)=",i15)') KZ_all,SUM(recvcnts)
      End If
      Allocate(bufrecv(KZ_all))
      Call mpi_allgatherv(bufsent, buffer_size, MPI_DOUBLE_COMPLEX, bufrecv, recvcnts, rdispls, MPI_DOUBLE_COMPLEX, tribeCOMM, mpi_err)

      ! Reconstruct densities one block at a time
      icount=0
      Do i=0,numberHFODDproc-1
         ! scalar densities
         Do IZ=KZ_min(i)+1,KZ_min(i)+block_vec(i)
            Do IY=1,NYHERM/IHABOY
               Do IX=1,NXHERM/IHABOX
                  icount=icount+1
                  DE_RHO(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                  icount=icount+1
                  DE_TAU(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                  icount=icount+1
                  DE_LPR(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                  icount=icount+1
                  DE_DIV(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                  icount=icount+1
                  PD_RHO(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                  icount=icount+1
                  PP_RHO(IX,IY,IZ,ITPNMX)=bufrecv(icount)
                End Do
            End Do
         End Do
         ! vector densities
         Do K=1,NDKART
            Do IZ=KZ_min(i)+1,KZ_min(i)+block_vec(i)
               Do IY=1,NYHERM/IHABOY
                  Do IX=1,NXHERM/IHABOX
                     icount=icount+1
                     DE_SPI(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_KIS(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_LPS(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_ROS(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_ROC(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_GRR(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                     icount=icount+1
                     DE_CUR(IX,IY,IZ,K,ITPNMX)=bufrecv(icount)
                  End Do
               End Do
            End Do
         End Do
         ! tensor densities
         Do L=1,NDKART
            Do K=1,NDKART
               Do IZ=KZ_min(i)+1,KZ_min(i)+block_vec(i)
                  Do IY=1,NYHERM/IHABOY
                     Do IX=1,NXHERM/IHABOX
                        icount=icount+1
                        DE_SCU(IX,IY,IZ,K,L,ITPNMX)=bufrecv(icount)
                        icount=icount+1
                        DE_DES(IX,IY,IZ,K,L,ITPNMX)=bufrecv(icount)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

#endif

      ! ATTENTION: In the versions between 2.22n and 2.37t, array "WARIGH" below was not multiplied by the
      !            occupation probabilities. As a result, the density matrix calculated by "DENMAC" was
      !            incorrect for IDEVAR=1. This bug was corrected in version 2.37u on 13/12/08.
      ! ATTENTION: In the versions 2.49l and x.xxx, array "WARIGH" was incorrectly multiplied by the
      !            occupation probabilities at finite-temperature.
      If(IDEVAR == 1 .And. IPAHFB < 1) Then
         Do ISPIN=0,NDSPIN
            Do ISTATE=1,LSTATE
               Do IBASE=1,LDBASE
                  WARIGH(IBASE,ISTATE,ISPIN)=WARIGH(IBASE,ISTATE,ISPIN)*V1OCCU(ISTATE)
               End Do
            End Do
         End Do
      End If ! End If IDEVAR

      !  Testing for the presence of any of the three plane symmetries. SYMDEN and SYMPAI can be called
      !  only when the calculations for the whole box are performed. This can be done by setting the
      !  parameter ISYMDE=1 in the input data.
      Allocate(DENALL(NDXHRM,NDYHRM,NDZHRM,NDALLD))
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,1) = DE_RHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,2) = DE_TAU(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,3) = DE_LPR(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,4) = DE_DIV(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM, 4           +1: 4+   NDKART)  = DE_SPI(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+   NDKART)+1:(4+ 2*NDKART)) = DE_KIS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 2*NDKART)+1:(4+ 3*NDKART)) = DE_GRR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 3*NDKART)+1:(4+ 4*NDKART)) = DE_LPS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 4*NDKART)+1:(4+ 5*NDKART)) = DE_ROS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 5*NDKART)+1:(4+ 6*NDKART)) = DE_ROC(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 6*NDKART)+1:(4+ 7*NDKART)) = DE_CUR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 7*NDKART)+1:(4+ 8*NDKART)) = DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 8*NDKART)+1:(4+ 9*NDKART)) = DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,2,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 9*NDKART)+1:(4+10*NDKART)) = DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,3,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+10*NDKART)+1:(4+11*NDKART)) = DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+11*NDKART)+1:(4+12*NDKART)) = DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,2,ITPNMX)
      DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+12*NDKART)+1:(4+13*NDKART)) = DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,3,ITPNMX)

      !If(PRINIT .And. ISYMDE == 1) Then
      !   Call SYMDEN(DENALL,NXHERM,NYHERM,NZHERM,NAMEPN)
      !   If(IPAHFB == 1) Call SYMPAI(ITPNMX,NXHERM,NYHERM,NZHERM)
      !End If

      ! The uncalculated parts of the box are here copied from the calculated parts. This is done according
      ! to the tables of symmetries defined explicitly in COPDEN and COPPAI.
      If(DOSYME) Then
         !Call COPDEN(DENALL,NXHERM,NYHERM,NZHERM,ISIGNY,ISIMPY,ISIQTY,ISIMTX,JSIMTY,ISIMTZ)
         !If(IPAHFB == 1) Call COPPAI(ITPNMX,NXHERM,NYHERM,NZHERM,IKERNE,ISIGNY,ISIMPY,ISIQTY,ISIMTX,JSIMTY,ISIMTZ)
         DE_RHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,1)
         DE_TAU(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,2)
         DE_LPR(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,3)
         DE_DIV(1:NDXHRM,1:NDYHRM,1:NDZHRM,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,4)
         DE_SPI(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM, 4          +1: 4+  NDKART)
         DE_KIS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+  NDKART)+1:(4+2*NDKART))
         DE_GRR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+2*NDKART)+1:(4+3*NDKART))
         DE_LPS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+3*NDKART)+1:(4+4*NDKART))
         DE_ROS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+4*NDKART)+1:(4+5*NDKART))
         DE_ROC(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+5*NDKART)+1:(4+6*NDKART))
         DE_CUR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+6*NDKART)+1:(4+7*NDKART))
         DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 7*NDKART)+1:(4+ 8*NDKART))
         DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,2,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 8*NDKART)+1:(4+ 9*NDKART))
         DE_SCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,3,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+ 9*NDKART)+1:(4+10*NDKART))
         DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+10*NDKART)+1:(4+11*NDKART))
         DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,2,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+11*NDKART)+1:(4+12*NDKART))
         DE_DES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,3,ITPNMX) = DENALL(1:NDXHRM,1:NDYHRM,1:NDZHRM,(4+12*NDKART)+1:(4+13*NDKART))
         Deallocate(DENALL)
      End If

      !Call CPUTIM('DENSHF',0)

   End Subroutine DENSHF_base
   !---------------------------------------------------------------------
   !> This subroutine simply prints an error message in case allocation of an array failed in a routine
   !---------------------------------------------------------------------
   Subroutine NOALLO(NAMMAT,NAMSUB)
      Character(Len=6), Intent(In) :: NAMMAT,NAMSUB

      Write(*,'(/,1X,19("/")," ALLOCATION OF MEMORY FOR ARRAY: ",A6,1X,19("/"),/, &
           &      1X,19("/")," WAS UNSUCCESSFUL IN SUBROUTINE: ",A6,1X,19("/"),/, &
           &      1X,19("/")," PLEASE RERUN THE CODE USING A COMPUTER",1X,19("/"),/, &
           &      1X,19("/")," THAT HAS A LARGER MEMORY. YOU MAY ALSO",1X,19("/"),/, &
           &      1X,19("/")," NEED TO USE A 64-BIT FORTRAN COMPILER.", 1X,19("/"),/)') NAMMAT,NAMSUB

      Stop ' ALLOCATION OF MEMORY FAILED'

   End Subroutine NOALLO
   !=======================================================================
   !
   !=======================================================================
End Module densities


!-----------------------------------------------------------------------
!
!                               MAIN PROGRAM
!
!-----------------------------------------------------------------------
Program test_denshf

   Use hfodd_types
   Use densities
   Use HObasis
#if(USE_MPI==1)
   Use MPI_structure
#endif

   Implicit None

   Logical :: PRINIT
#if(USE_MPI==1)
   Integer(ipr) :: mpi_err
#endif
   Integer(ipr) :: j,KX,KY,KZ,speed
   Integer(ipr) :: NGAUSS,NOSCIL,NLIMIT,IPAHFB,MREVER,ICHARG,KARTEZ,NOBODY
   Integer(ipr) :: ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY
   Integer(ipr) :: IPNMIX,ITPNMX,ISYMDE,IDEVAR,IKERNE
   Real(pr) :: ELIMIT,HOSCAX,HOSCAY,HOSCAZ,X_MASS,DENTOT,W_HERM,cput_start
   Character(Len=8) :: NAMEPN

   ! If MPI is used, define the subcommunicators
#if(USE_MPI==1)
   Call create_MPI_layout(cput_start)
#endif

   ! Allocate all arrays needed in the following
   Call allocate_data()

   !=======================================================================
   !         CHOOSE SPEED OF THE CODE BELOW (1 - SLOW to 3 - FAST)
   speed = 1
   Select Case (speed)
   Case (1)
      ! Slow (~ 2 min with serial, 4 threads)
      NOSCIL=30; NLIMIT=2500; NXHERM=46; NYHERM=46; NZHERM=66
      NOSCIL=60; !NLIMIT=5000; NXHERM=92; NYHERM=92; NZHERM=132
   Case (2)
      ! Not so fast (~ 20 sec with serial, 4 threads)
      NOSCIL=18; NLIMIT=1330; NXHERM=40; NYHERM=40; NZHERM=40
   Case (3)
      ! Very fast (~ 1 sec with serial, 4 threads)
      NOSCIL=8; NLIMIT=165; NXHERM=20; NYHERM=20; NZHERM=20
   Case Default
      NOSCIL=8; NLIMIT=165; NXHERM=20; NYHERM=20; NZHERM=20
   End Select
   !=======================================================================

   NGAUSS=NZHERM

   ! Define HO basis functions: arrays CERMTS(:,:,:), CHRMTS(:,:,:) and CDHRMT(:,:,:)
   NOBODY=2; HOSCAX=0.5_pr; HOSCAY=0.5_pr; HOSCAZ=0.35_pr
   Call DEFINT(NXHERM,NYHERM,NZHERM,NOBODY,HOSCAX,HOSCAY,HOSCAZ)
   KARTEZ=1
   Call COPHER(KARTEZ,NGAUSS)
   KARTEZ=2
   Call COPHER(KARTEZ,NGAUSS)
   KARTEZ=3
   Call COPHER(KARTEZ,NGAUSS)

   ! Define basis states: arrays INDICE(:,:,:), LAYOFX(:), etc., overall size LDBASE
   ELIMIT=0.0_pr
   Call SETBAS(NOSCIL,NLIMIT,ELIMIT)

   ! Preset quantities (do not change)
   LDTOTS = LDBASE
   LDSTAT = LDBASE
   LDUPPE = LDBASE
   LDTIMU = LDBASE
   Do j=1,LDBASE
      WARIGH(j,j,0)=Sqrt(0.05)
      WALEFT(j,j,0)=Sqrt(0.95)
   End Do

   ! Deactivating symmetrization (= not arbitrarily reducing size of loops)
   IHABOX = 1
   IHABOY = 1
   IHABOZ = 1
   DOSYME = .False.

   ! Calling the routine
   ISIMTX=1
   JSIMTY=1
   ISIMTZ=1
   ISIGNY=1
   ISIMPY=1
   ISIQTY=1
   IPAHFB=1
   MREVER=0
   ICHARG=0
   IPNMIX=0
   ITPNMX=ICHARG
   IDEVAR=1
   ISYMDE=0
   IKERNE=0

   ! Calling the routine
   PRINIT = .True.
   NAMEPN = 'xxxxxxxx'
   Call DENSHF_base(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ, &
                    ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,ICHARG, &
                    IPNMIX,ITPNMX, &
                    NXMAXX,NAMEPN,PRINIT,IDEVAR,ISYMDE,IKERNE)

   ! Verifying the result: whatever modifications you make, you must always get these numbers at the end
#if(USE_MPI==1)
   If(worldRank == 0) Then
#endif
      DENTOT=0.0_pr; X_MASS=0.0_pr
      Do KZ=1,NZHERM
         Do KY=1,NYHERM
            Do KX=1,NXHERM
               DENTOT=REAL(DE_RHO(KX,KY,KZ,0))*EXPAUX(KX,1)*EXPAUX(KY,2)*EXPAUX(KZ,3)
               W_HERM=WGTACT(KX,1)*WGTACT(KY,2)*WGTACT(KZ,3)
               X_MASS=X_MASS+W_HERM*DENTOT
            End Do
         End Do
         Write(*,'(I4,2E24.10)') KZ,DE_RHO(1,1,KZ,0)*EXPAUX(1,1)*EXPAUX(1,2)*EXPAUX(KZ,3)
      End Do
      Write(*,'(E24.10)') X_MASS
#if(USE_MPI==1)
   End If
#endif

   ! Deallocating arrays properly
   Call deallocate_data()

#if(USE_MPI==1)
   ! Close MPI environment
   Call mpi_finalize(mpi_err)
#endif

End Program test_denshf
