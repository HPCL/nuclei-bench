!***********************************************************************
!
!    Copyright (c) 2020, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Lead developer: Nicolas Schunck, schunck1@llnl.gov
!    HFODD
!    -----
!      LLNL-CODE-710577 All rights reserved.
!      LLNL-CODE-470611 All rights reserved.
!
!      Copyright 2017, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
!                      J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
!                      X.B. Wang, and T.R. Werner
!      Copyright 2012, N. Schunck, J. Dobaczewski, J. McDonnell,
!                      W. Satula, J.A. Sheikh, A. Staszczak,
!                      M. Stoitsov, P. Toivanen
!      Copyright 2009, J. Dobaczewski, W. Satula, B.G. Carlsson, J. Engel,
!                      P. Olbratowski, P. Powalowski, M. Sadziak,
!                      J. Sarich, N. Schunck, A. Staszczak, M. Stoitsov,
!                      M. Zalewski, H. Zdunczuk
!      Copyright 2004, 2005, J. Dobaczewski, P. Olbratowski
!      Copyright 1997, 2000, J. Dobaczewski, J. Dudek
!
!    HFBTHO
!    -----
!      LLNL-CODE-728299 All rights reserved.
!      LLNL-CODE-573953 All rights reserved.
!
!      Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang, J. Sarich
!      Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                      N. Michel, J. Sarich, S. Wild
!      Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!
!    This file is part of DFTNESS.
!
!    DFTNESS is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    DFTNESS is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with DFTNESS. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

 !-------------------------------------------------------------------
 !
 !                  TYPES, CONSTANTS
 !
 !-------------------------------------------------------------------
Module hfodd_types

   Implicit None

   Integer, Parameter, Public :: ipr=Kind(1)     !< default integer precision
   Integer, Parameter, Public :: pr =Kind(1.0d0) !< default real precision

#if(USE_MPI==0)
   Integer, Parameter, Public :: NFIPRI = 6
#else
   Integer, Parameter, Public :: NFIPRI = 88
#endif

   Real(pr), Parameter, Public :: PI = 4.0_pr*Atan(1.0_pr) !< Pi
   Real(pr), Parameter, Public :: DEGRAD = 180.0_pr/PI     !< To convert from radians to degrees
   Real(pr), Parameter, Public :: SQR4PI = 4.0_pr*Sqrt(Atan(1.0_pr)) !< \f$ \sqrt{4\pi} \f$

   Complex(pr), Parameter, Public :: C_ZERO=Cmplx( 0.0_pr,0.0_pr)
   Complex(pr), Parameter, Public :: C_UNIT=Cmplx( 1.0_pr,0.0_pr)
   Complex(pr), Parameter, Public :: UNIT_I=Cmplx( 0.0_pr,1.0_pr)
   Complex(pr), Parameter, Public :: C_MINU=Cmplx(-1.0_pr,0.0_pr)

   ! Physical constants
   Real(pr), Parameter, Public :: H_BARC=197.32891_pr, XMASSP=938.27231_pr, XMASSN=938.90590_pr
   Real(pr), Parameter, Public :: HBCOE2=137.03602_pr

   ! (Possibly) convenient derived types
   Type nucleus
      Integer(ipr), Public :: IZ_FIX
      Integer(ipr), Public :: IN_FIX
      Real(pr), Public :: XZ_FIX
      Real(pr), Public :: XN_FIX
   End Type nucleus

   Type skyrme_contributions
      Real(pr) :: ERHO_T
      Real(pr) :: ERHO_S
      Real(pr) :: ERHODT
      Real(pr) :: ERHODS
      Real(pr) :: ELPR_T
      Real(pr) :: ELPR_S
      Real(pr) :: ETAU_T
      Real(pr) :: ETAU_S
      Real(pr) :: ESCU_T
      Real(pr) :: ESCU_S
      Real(pr) :: EDIV_T
      Real(pr) :: EDIV_S
      Real(pr) :: ESPI_T
      Real(pr) :: ESPI_S
      Real(pr) :: ESPIDT
      Real(pr) :: ESPIDS
      Real(pr) :: ELPS_T
      Real(pr) :: ELPS_S
      Real(pr) :: ECUR_T
      Real(pr) :: ECUR_S
      Real(pr) :: EKIS_T
      Real(pr) :: EKIS_S
      Real(pr) :: EROT_T
      Real(pr) :: EROT_S
   End Type skyrme_contributions

   Type energies
      Real(pr) :: EKIN_N
      Real(pr) :: EKIN_P
      Real(pr) :: EKIN_T
      Real(pr) :: EPOT_N
      Real(pr) :: EPOT_P
      Real(pr) :: EPOT_T
      Real(pr) :: ESUM_N
      Real(pr) :: ESUM_P
      Real(pr) :: ESUM_T
      Real(pr) :: EPAI_N
      Real(pr) :: EPAI_P
      Real(pr) :: EPAI_T
      Real(pr) :: EREA_N
      Real(pr) :: EREA_P
      Real(pr) :: EREA_T
      Real(pr) :: ELIP_N
      Real(pr) :: ELIP_P
      Real(pr) :: ELIP_T
      Real(pr) :: ECOULD
      Real(pr) :: ECOULE
      Real(pr) :: ECOULT
      Real(pr) :: ECOULS
      Real(pr) :: ECOULV
      Real(pr) :: EMULCO
      Real(pr) :: EMUSLO
      Real(pr) :: EMUREA
      Real(pr) :: ESIFCO
      Real(pr) :: ESISLO
      Real(pr) :: ESIREA
      Real(pr) :: ESPICO
      Real(pr) :: ESPSLO
      Real(pr) :: ESPREA
      Real(pr) :: ENREAR
      Real(pr) :: ECORCM
      Real(pr) :: ECOR_R
      Real(pr) :: EEVEW0
      Real(pr) :: EODDW0
      Real(pr) :: ENE_W0
      Real(pr) :: ENEVEN
      Real(pr) :: ENEODD
      Real(pr) :: ENESKY
      Real(pr) :: ESTABN
      Real(pr) :: ETOTSP
      Real(pr) :: ETOTFU
   End Type energies

Contains

End Module hfodd_types