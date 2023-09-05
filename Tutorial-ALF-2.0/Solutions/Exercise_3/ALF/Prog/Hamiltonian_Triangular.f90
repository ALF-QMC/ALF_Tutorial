!  Copyright (C) 2016 - 2018 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de 
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing   
!       to the ALF project or to mark your material in a reasonable way as different from the original version 


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> This module defines the  Hamiltonian and observables.  Here, we have included a
!> set of predefined Hamiltonians. They include the Hubbard and SU(N) tV models
!> on honeycomb, pi-flux and square lattices.

!> @details
!> The public variables of this module are the following
!>
!> 
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable 
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!> 
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable  
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L   
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma(:,:) 
!> \verbatim Type(Fields)
!> Array containing all auxiliary fields. The first index runs through the operator sequence. The second
!> through the time slies.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!> 
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot  
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!> 
!> @param [public] Group_Comm 
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines. 
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>  
!> You still have to add some docu for the other private variables in this module.      
!>
!--------------------------------------------------------------------

    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Predefined_structures
      Use Fields_mod
      Use LRC_Mod

      
      Implicit none

     
      Type (Operator),     dimension(:,:), allocatable :: Op_V 
      Type (Operator),     dimension(:,:), allocatable :: Op_T
      Type (WaveFunction), dimension(:),   allocatable :: WF_L
      Type (WaveFunction), dimension(:),   allocatable :: WF_R
      Type  (Fields)       :: nsigma
      Integer              :: Ndim
      Integer              :: N_FL
      Integer              :: N_SUN
      Integer              :: Ltrot
      Integer              :: Thtrot 
      Logical              :: Projector
      Integer              :: Group_Comm
      Logical              :: Symm


      Type (Lattice),       private :: Latt
      Type (Unit_cell),     private :: Latt_unit
      Integer,              private :: L1, L2, N_curr_bonds, ham_np
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_tV, ham_t2, ham_t3, ham_t4
      real (Kind=Kind(0.d0)),        private :: ham_lambda, Ham_h, Ham_J, Ham_xi, ham_cdw, ham_b, ham_v
      real (Kind=Kind(0.d0)),        private :: ham_alpha, Percent_change
      real (Kind=Kind(0.d0)),        private :: XB_Y, Phi_Y, XB_X, Phi_X
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional, Checkerboard, ham_flux, chem_tune, Mz_decoupling
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),  private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer,  allocatable ,  private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)
      Real (Kind=Kind(0.d0)),  private :: Bond_pos(2,12)

      
    contains 

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set
#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: ierr, N_part, nf, nc, n, I, no, I1, nc1, J1, I2, J2, I_tmp
          Character (len=64) :: file_info, file_para
          Complex(kind=kind(0.d0)) :: flux
          Real (kind=kind(0.d0))   :: x_p(2), signum, s1, s2, BC, factor
          COMPLEX(Kind=Kind(0.d0)), allocatable :: H0(:,:), U(:,:)
          Real(Kind=Kind(0.d0)), allocatable :: En(:)

          
          
          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter. 
          
          NAMELIST /VAR_Lattice/  L1, L2, Lattice_type, Model,  Checkerboard, N_SUN, Phi_X, XB_X, Symm

          NAMELIST /VAR_Triangular/  ham_T, ham_chem, ham_np, Dtau, Beta, Theta, Projector, ham_V

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
#endif
          ! Global "Default" values.
          N_SUN        = 1
          Checkerboard = .false.
          Symm         = .false.
          Phi_X        = 0.d0
          XB_X         = 1.d0
          Phi_Y        = 0.d0
          XB_Y         = 1.d0

#ifdef MPI
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          !if ( irank_g == 0 )   write(6,*) "Mpi Test", igroup, isize_g
#endif
          ! Open files
#if defined(MPI) 
          If (Irank_g == 0 ) then
#endif
             File_para = "parameters"
             File_info = "info"
#if defined(TEMPERING) 
             write(File_para,'(A,I0,A)') "Temp_",igroup,"/parameters"
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             OPEN(Unit = 50,file=file_info,status="unknown",position="append")
#ifdef MPI
          Endif
#endif


#ifdef MPI
          If (Irank_g == 0 ) then
#endif
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
 
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Phi_X       ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(XB_X        ,1  ,MPI_REAL8  ,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
          CALL MPI_BCAST(Checkerboard,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Symm        ,1  ,MPI_LOGICAL  , 0,Group_Comm,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,Group_Comm,IERR)
#endif
          
          Call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)
          
#if defined(MPI) 
          If (Irank == 0 ) then
#endif
            OPEN(UNIT=123,FILE='Triangular_Lattice',STATUS='replace',ACTION='write',IOSTAT=ierr)
            Do I=1,Ndim
              x_p = Latt%list(List(I,1),1)*Latt%a1_p + Latt%list(List(I,1),2)*Latt%a2_p
              x_p = x_p + Latt_Unit%Orb_pos_p(List(I,2),:)
              write(123,*) I, List(I,1), List(I,2),  Latt%list(List(I,1),1), Latt%list(List(I,1),2), x_p(1), x_p(2)
            enddo
#if defined(MPI)
          endif
#endif


#ifdef MPI
          If (Irank_g == 0) then
#endif
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) '# of orbitals : ', Ndim
             Write(50,*) 'Checkerboard  : ', Checkerboard
             Write(50,*) 'Symm. decomp  : ', Symm
#ifdef MPI
          Endif
#endif



          ! Default is finite temperature. 
          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
          N_FL  = 1
          N_SUN = 1
          ham_t = 1.0d0
          ham_V = 0.0d0
          ham_flux = .false.
          ham_chem = 0.d0
             ham_np=Ndim*N_FL*N_SUN/2
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                READ(5,NML=VAR_Triangular)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                  Write(50,*) 'N_part        : ', ham_np
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_V         : ', Ham_V
                Write(50,*) 'Ham_chem      : ', Ham_chem
                if (ham_flux) Write(50,*) 'Pi-flux present!'
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot      ,1,MPI_INTEGER,0,Group_Comm,ierr)
             CALL MPI_BCAST(Thtrot     ,1,MPI_INTEGER,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_np     ,1,MPI_INTEGER,0,Group_Comm,ierr)
             CALL MPI_BCAST(Projector  ,1,MPI_LOGICAL,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_flux   ,1,MPI_LOGICAL,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_T      ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(ham_chem   ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(Dtau       ,1,MPI_REAL8  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(Beta       ,1,MPI_REAL8  ,0,Group_Comm,ierr)
#endif
            N_curr_bonds=6
            if (ham_t3 .ne. 0.0d0) N_curr_bonds=10
            if (ham_t4 .ne. 0.0d0) N_curr_bonds=12
            
            Bond_pos      = 0.00d0
            
            Bond_pos(2,2) = 0.50d0
            Bond_pos(2,3) =-0.50d0
            Bond_pos(2,4) = 1.00d0
            Bond_pos(2,5) = 1.50d0
            Bond_pos(2,6) = 0.50d0
                  
            Bond_pos(1,7) = 0.50d0
            Bond_pos(2,7) = 1.00d0
            Bond_pos(1,8) = 0.50d0
            Bond_pos(2,8) =-1.00d0
            Bond_pos(1,9) = 0.50d0
            Bond_pos(2,9) = 2.00d0
            Bond_pos(1,10)= 0.50d0
            Bond_pos(2,10)= 0.00d0
            
            Bond_pos(1,11)= 0.50d0
            Bond_pos(2,11)= 0.00d0
            Bond_pos(1,12)= 0.50d0
            Bond_pos(2,12)= 1.00d0
            
            allocate(Op_T(1,N_FL))
            do n = 1,N_FL
              Call Op_make(Op_T(1,n),Ndim)
              flux = 1.0d0
              if (ham_flux) flux = cmplx(1.d0, 1.d0, kind(0.d0))/sqrt(2.d0)
              if (n==2) flux = conjg(flux)
              nc = 1
              DO I = 1, Latt%N
                  do no = 1,Latt_unit%Norb
                    I1 = Invlist(I,no)
                    if (no==1) then
                      Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    else
                      Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    endif
                  enddo
                  I1 = Invlist(I,1)
                  J1 = I1
                  Do nc1 = 1,Latt_unit%N_coord
                    select case (nc1)
                    case (1)
                        J1 = invlist(Latt%nnlist(I,0, 1),1) 
                    case (2)
                        J1 = invlist(Latt%nnlist(I,1, 0),1) 
                    case (3)
                        J1 = invlist(Latt%nnlist(I,1,-1),1) 
                    case default
                        Write(6,*) ' Error in  Ham_Hop '  
                        Stop
                    end select
                    write(*,*) I1, J1
                    Op_T(nc,n)%O(I1,J1) = Op_T(nc,n)%O(I1,J1) - Ham_T
                    Op_T(nc,n)%O(J1,I1) = Op_T(nc,n)%O(J1,I1) - Ham_T
                  Enddo
              Enddo
              Do I = 1,Ndim
                Op_T(nc,n)%P(i) = i 
              Enddo
              if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6) then 
                Op_T(nc,n)%g = 0.d0
              else
                Op_T(nc,n)%g = -Dtau
              endif
              Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
              Call Op_set(Op_T(nc,n))
!               write(*,*)
!               Do I = 1,Size(Op_T(nc,n)%E,1)
!                 Do J1 = 1,Size(Op_T(nc,n)%E,1)
!                   if (abs(Op_T(nc,n)%O(I,J1)) > 1.E-6 ) write(*,*) I, J1, real(Op_T(nc,n)%O(I,J1)), aimag(Op_T(nc,n)%O(I,J1))
!                 Enddo
!               Enddo
!               write(*,*)
!               Write(6,*) Op_T(nc,n)%E(:)
!               write(*,*)
            Enddo
            
#ifdef MPI
            If (Irank_g == 0) then
#endif
              write(50,*) "Gap is:              ", Op_t(1,1)%E(Ndim/2+1) - Op_t(1,1)%E(Ndim/2)
              write(50,*) "Bandwidth is:        ", Op_t(1,1)%E(Ndim/2) - Op_t(1,1)%E(1)
              write(50,*) "ratio width/gap is:  ", (Op_t(1,1)%E(Ndim/2)-Op_t(1,1)%E(1))/&
                                                  &(Op_t(1,1)%E(Ndim/2+1)-Op_t(1,1)%E(Ndim/2))
!               write(*,*) Op_t(1,1)%E(Ndim/4),Op_t(1,1)%E(Ndim/4+1)
#ifdef MPI
            Endif
#endif             
          
          
          if (Projector) then
             Allocate(WF_L(N_FL),WF_R(N_FL))
             N_part = ham_np/N_FL/N_SUN
             factor=1.d0
             if (N_part==2) factor=0.999d0
             BC=1.d0
             if (mod(L1,8)==0 .and. mod(L2,8)==0 .and. N_part.ne.2) BC=-1.d0
             Allocate(H0(Ndim,Ndim),U(Ndim,Ndim),En(Ndim))
             Do nf = 1,N_FL
                Call WF_alloc(WF_L(nf),Ndim,N_part)
                Call WF_alloc(WF_R(nf),Ndim,N_part)
                
                if (N_part==1) then
                  WF_L(nf)%P(:,1)=1.d0/sqrt(dble(Ndim))
                  WF_R(nf)%P(:,1)=1.d0/sqrt(dble(Ndim))
                  WF_L(nf)%Degen = 1.d0
                  WF_R(nf)%Degen = 1.d0
                else
                  H0=0.d0
                  
                  flux = 1.0d0
                  if (ham_flux) flux = cmplx(1.d0, 1.d0, kind(0.d0))/sqrt(2.d0)
                  if (nf==2) flux = conjg(flux)
                  nc = 1
                  DO I = 1, Latt%N
                    do no = 1,Latt_unit%Norb
                      I1 = Invlist(I,no)
                      H0(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    Do nc1 = 1,Latt_unit%N_coord
                      select case (nc1)
                      case (1)
                          J1 = invlist(I,2) 
                      case (2)
                          J1 = invlist(Latt%nnlist(I,0, 1),2) 
                      case (3)
                          J1 = invlist(Latt%nnlist(I,-1,1),2) 
                      case (4)
                          J1 = invlist(Latt%nnlist(I,-1,0),2) 
                      case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                      end select
                      signum=1.d0
                      ! x-coordinate is zero
                      if ( Latt%list(I,1)+Latt%list(I,2)==0 .and. nc1==2 ) signum=BC
                      ! x-coordinate is one
                      if ( Latt%list(I,1)+Latt%list(I,2)==1 .and. nc1==4 ) signum=BC
                      ! y-coordinate is zero
                      if ( Latt%list(I,1)-Latt%list(I,2)==0 .and. nc1==1 ) signum=BC
                      ! y-coordinate is one
                      if ( Latt%list(I,1)-Latt%list(I,2)==1 .and. nc1==3 ) signum=BC
                      if (nc1 == 1 .or. nc1 == 3 ) then
                          ! odd bond
                          H0(I1,J1) = H0(I1,J1) - signum*Ham_T*flux
                          H0(J1,I1) = H0(J1,I1) - signum*Ham_T*conjg(flux)
                      Else
                          ! even bond
                          H0(I1,J1) = H0(I1,J1) - signum*Ham_T*conjg(flux)
                          H0(J1,I1) = H0(J1,I1) - signum*Ham_T*flux
                      endif
                    Enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    I2 = Invlist(I,2)
                    J2 = I2
                    Do nc1 = 1,2
                      select case (nc1)
                      case (1)
                          J1 = invlist(Latt%nnlist(I,1, 0),1)
                          J2 = invlist(Latt%nnlist(I,1, 0),2) 
                      case (2)
                          J1 = invlist(Latt%nnlist(I,0, 1),1)
                          J2 = invlist(Latt%nnlist(I,0, 1),2)
                      case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                      end select
                      signum=1.d0
                      s2=1.d0
                      s1=1.d0
                      ! x-coordinate is zero
                      if ( Latt%list(I,1)+Latt%list(I,2)==0 ) signum=BC
                      if ( Latt%list(I,1)-Latt%list(I,2)==0 .and. nc1==1 ) s1=BC
                      if ( Latt%list(I,1)-Latt%list(I,2)==1 .and. nc1==2 ) s1=BC
                      if ( Latt%list(I,1)-Latt%list(I,2)==0 .and. nc1==2 ) s2=BC
                      if ( Latt%list(I,1)-Latt%list(I,2)==-1.and. nc1==1 ) s2=BC
                      if (nc1 == 1) then
                          ! odd bond
                          H0(I1,J1) = H0(I1,J1) - signum*s1*Ham_T2*factor
                          H0(J1,I1) = H0(J1,I1) - signum*s1*Ham_T2*factor
                          H0(I2,J2) = H0(I2,J2) + signum*s2*Ham_T2*factor
                          H0(J2,I2) = H0(J2,I2) + signum*s2*Ham_T2*factor
                      Else
                          ! even bond
                          H0(I1,J1) = H0(I1,J1) + signum*s1*Ham_T2*factor
                          H0(J1,I1) = H0(J1,I1) + signum*s1*Ham_T2*factor
                          H0(I2,J2) = H0(I2,J2) - signum*s2*Ham_T2*factor
                          H0(J2,I2) = H0(J2,I2) - signum*s2*Ham_T2*factor
                      endif
                    Enddo
                  enddo
                  
                  Call Diag(H0,U,En)

                  do I2=1,N_part
                    do I1=1,Ndim
                        WF_L(nf)%P(I1,I2)=U(I1,I2)
                        WF_R(nf)%P(I1,I2)=U(I1,I2)
                    enddo
                  enddo
                  WF_L(nf)%Degen = En(N_part+1) - En(N_part)
                  WF_R(nf)%Degen = En(N_part+1) - En(N_part)
                endif
             enddo
             deallocate(H0,U,En)

#ifdef MPI
             If (Irank_g == 0) then
#endif
                Do nf = 1,N_FL
                   Write(50,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(50,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
                   
#ifdef MPI
             Endif
#endif             
             
          endif

#ifdef MPI
          If (Irank_g == 0 )  then
#endif
             close(50)
             Close(5)
#ifdef MPI
          endif
#endif
          
          call Ham_Vint
          


        end Subroutine Ham_Set
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
        Subroutine Ham_Vint
          
          Implicit none 
          
          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X

          Allocate(Op_V(Latt_unit%N_coord*Ndim,N_FL))
          do nf = 1,N_FL
            do i  = 1, Ndim
                Call Op_make(Op_V(i,nf),2)
            enddo
          enddo
          
          Do nf = 1,N_FL
            nc = 0
              Do i = 1,Latt%N
                I1 = Invlist(I,1)
                I2 = I1
                
                Do nc1 = 1,Latt_unit%N_coord
                  select case (nc1)
                  case (1)
                      I2 = invlist(Latt%nnlist(I,0, 1),1) 
                  case (2)
                      I2 = invlist(Latt%nnlist(I,1, 0),1) 
                  case (3)
                      I2 = invlist(Latt%nnlist(I,1,-1),1) 
                  case default
                      Write(6,*) ' Error in  Ham_Hop '  
                      Stop
                  end select
                    
                  nc = nc + 1
                  Call Op_make(Op_V(nc,nf),2)
                  Op_V(nc,nf)%P(1) = I1
                  Op_V(nc,nf)%P(2) = I2
                  
                  Op_V(nc,nf)%alpha  = cmplx(0.0d0, 0.d0, kind(0.D0))
                  select case (model)
                  case ("den")
                    Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(2,2) = cmplx(1.d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%alpha  = cmplx(-1.0d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_V/2.0d0, 0.D0, kind(0.D0))) 
                  case ("mz")
                    Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(2,2) = cmplx(-1.d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%g      = SQRT(CMPLX(+DTAU*ham_V/2.0d0, 0.D0, kind(0.D0))) 
                  case ("hop")
                    Op_V(nc,nf)%O(1,2) = cmplx(1.d0, 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(2,1) = cmplx(1.d0, 0.d0, kind(0.D0)) 
                    Op_V(nc,nf)%g      = SQRT(CMPLX(+DTAU*ham_V/2.0d0, 0.D0, kind(0.D0))) 
                  case default
                    Write(6,*) ' Error in  Ham_V '  
                    Stop
                  end select
                  Op_V(nc,nf)%type   = 2
                  Call Op_set( Op_V(nc,nf) )
                Enddo
              Enddo
          Enddo
          
        end Subroutine Ham_Vint

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
        Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)  
          Implicit none
          !> Operator index
          Integer, Intent(IN) :: n
          !> Time slice
          Integer, Intent(IN) :: nt
          !> New local field on time slice nt and operator index n
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new
          
          Integer :: nt1,I
          !Write(6,*) "Hi1"
          
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then 
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma%i(n,nt)*nsigma%i(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1 
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             nt1 = nt - 1 
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif

          If (model == "LRC" ) then
             !Write(6,*) "Hi2"
             S0 = LRC_S0(n,dtau,nsigma%f(:,nt),Hs_new,N_SUN)
          Endif
          
        end function S0
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Setup  bonds and median lattice. 
!> @details
!> Note that the median lattice on which the Ising bond fields are defined could be defined generally in
!> the predefined lattice subroutine.
!--------------------------------------------------------------------

        Subroutine Setup_Ising_action
          
          ! This subroutine sets up lists and arrays so as to enable an 
          ! an efficient calculation of  S0(n,nt) 

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*Latt_unit%N_coord,2) )
          
          nc = 0
          do nth = 1,2*Latt_unit%N_coord  
             Do n1= 1, L1/2
                Do n2 = 1,L2
                   nc = nc + 1
                   n_orientation = 1
                   I1 = 1
                   If (nth == 1 ) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 2) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 3) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   elseif (nth == 4) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p  + latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   endif
                   L_bond(I1,n_orientation) = nc
                   L_bond_inv(nc,1) = I1  
                   L_bond_inv(nc,2) = n_orientation 
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                Enddo
             Enddo
          Enddo
          ! Setup the nearest neigbour lists for the Ising spins. 
          allocate(Ising_nnlist(2*Latt%N,4)) 
          do I  = 1,Latt%N
             n  = L_bond(I,1)
             n1 = L_bond(Latt%nnlist(I, 1, 0),2)
             n2 = L_bond(Latt%nnlist(I, 0, 0),2)
             n3 = L_bond(Latt%nnlist(I, 0,-1),2)
             n4 = L_bond(Latt%nnlist(I, 1,-1),2)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
             n  = L_bond(I,2)
             n1 = L_bond(Latt%nnlist(I, 0, 1),1)
             n2 = L_bond(Latt%nnlist(I,-1, 1),1)
             n3 = L_bond(Latt%nnlist(I,-1, 0),1)
             n4 = L_bond(Latt%nnlist(I, 0, 0),1)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
          enddo
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J) 
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J) 
          
        End Subroutine Setup_Ising_action
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No, Norb
          Character (len=64) ::  Filename


          Norb = Latt_unit%Norb
          ! Scalar observables
          no=1
!           if( projector ) no=0
          Allocate ( Obs_scal(6+no) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Kin"
             case (2)
                N = 1;   Filename ="Pot"
             case (3)
                N = 1;   Filename ="Part"
             case (4)
                N = 1;   Filename ="Ener"
             case (5)
                N = 1;   Filename ="Kin_x"
             case (6)
                N = 2;   Filename ="Doub_Occ"
             ! this has to remain in the last position
             case (7)
                N = 11;   Filename ="Jx_sW"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(14) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = Norb;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
             case (3)
                Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
             case (4)
                Ns = Latt%N;  No = Norb;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = Norb;  Filename ="SpinT"
             case (6)
                Ns = Latt%N;  No = Norb;  Filename ="SC"
             case (7)
                Ns = Latt%N;  No = N_curr_bonds;  Filename ="Curr"
             case (8)
                Ns = Latt%N;  No = 4;  Filename ="BondDen"
             case (9)
                Ns = Latt%N;  No = Norb;  Filename ="Green_up"
             case (10)
                Ns = Latt%N;  No = Norb;  Filename ="Green_down"
             case (11)
                Ns = Latt%N;  No = Norb;  Filename ="SC_sym2"
             ! This should remain the last entry
             case (12)
                Ns = Latt%N;  No = Norb;  Filename ="Den_sus"
             case (13)
                Ns = Latt%N;  No = Norb;  Filename ="Spin_sus"
             case (14)
                Ns = Latt%N;  No = N_curr_bonds;  Filename ="Curr_sus"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             No=0
             If (projector) no=1
             Allocate ( Obs_tau(8+no) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
!                 case (3)
!                    Ns = Latt%N; No = Norb;  Filename ="SpinXY"
!                 case (4)
!                    Ns = Latt%N; No = Norb;  Filename ="Den"
                case (2)!5
                   Ns = Latt%N; No = Norb;  Filename ="SC"
                case (3)!6
                   Ns = Latt%N; No = N_curr_bonds;  Filename ="Curr"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="Green_up"
                case (5)
                   Ns = Latt%N; No = Norb;  Filename ="Green_down"
                case (6)!5
                   Ns = Latt%N; No = Norb;  Filename ="SC_sym2"
                case (7)!5
                   Ns = Latt%N; No = Norb;  Filename ="Den"
                case (8)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                ! theo following case should remain the last one
                case (9)!7
                   Ns = Latt%N; No = Norb;  Filename ="Green_0T"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1-2*Thtrot
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        End Subroutine Alloc_obs
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Global moves
!> 
!> @details
!>  This routine generates a 
!>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!> @param [OUT]  T0_Proposal_ratio Real
!> \verbatimam
!>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)  
!> \endverbatim
!> @param [OUT]  Size_clust Real
!> \verbatim
!>  Size of cluster that will be flipped.
!> \endverbatim
!-------------------------------------------------------------------
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Type (Fields),  Intent(IN)  :: nsigma_old

          ! Local
          Integer :: N_op, N_tau, n1,n2, n

          T0_Proposal_ratio  = 1.d0
          size_clust         = 3.d0

          nsigma%f = nsigma_old%f
          nsigma%t = nsigma_old%t
          N_op  = size(nsigma_old%f,1)
          N_tau = size(nsigma_old%f,2)
          Do n  = 1, Nint(size_clust)
             n1 = nranf(N_op)
             n2 = nranf(N_tau)
             nsigma%f(n1,n2) = nsigma_old%flip(n1,n2)
          enddo
          
        End Subroutine Global_move
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Computes the ratio exp(S0(new))/exp(S0(old))
!> 
!> @details
!> This function computes the ratio \verbatim  e^{-S0(nsigma)}/e^{-S0(nsigma_old)} \endverbatim
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!-------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          ! Arguments
          Type (Fields),  INTENT(IN) :: nsigma_old

          ! Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m
         

          Delta_S0_global = 1.d0
          If ( Model == "Hubbard_SU2_Ising" ) then
             nc_F = 0
             nc_J = 0
             nc_h_p = 0
             nc_h_m = 0
             Do I = 1,Latt%N
                n1  = L_bond(I,1)
                n2  = L_bond(Latt%nnlist(I,1,0),2)
                n3  = L_bond(Latt%nnlist(I,0,1),1)
                n4  = L_bond(I,2)
                do nt = 1,Ltrot
                   nt1 = nt +1 
                   if (nt == Ltrot) nt1 = 1
                   if (nsigma%i(n1,nt) == nsigma%i(n1,nt1) ) then 
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n1,nt) == nsigma_old%i(n1,nt1) ) then 
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   if (nsigma%i(n4,nt) == nsigma%i(n4,nt1) ) then 
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n4,nt) == nsigma_old%i(n4,nt1) ) then 
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif
                   
                   nc_F = nc_F + nsigma%i    (n1,nt)*nsigma%i    (n2,nt)*nsigma%i    (n3,nt)*nsigma%i  (n4,nt)  &
                        &      - nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt) 
                   
                   nc_J = nc_J + nsigma%i(n1,nt)*nsigma%i(n2,nt) + &
                        &        nsigma%i(n2,nt)*nsigma%i(n3,nt) + &
                        &        nsigma%i(n3,nt)*nsigma%i(n4,nt) + &
                        &        nsigma%i(n4,nt)*nsigma%i(n1,nt) - &
                        &        nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt) - &
                        &        nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt) - &
                        &        nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt) - &
                        &        nsigma_old%i(n4,nt)*nsigma_old%i(n1,nt) 
                   
                enddo
             enddo
             !             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
             !                  &            exp( -Dtau*(Ham_F*real(nc_F,kind(0.d0)) -  Ham_J*real(nc_J,kind(0.d0))))
             ! No flux in example code. May want to include it.
             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
                  &            exp( Dtau* Ham_J*real(nc_J,kind(0.d0)))
             !Write(6,*) Delta_S0_global
          endif


        end Function Delta_S0_global
        
        
!==========================================================   

        Complex (Kind=kind(0.d0)) Function   Corr(GR, GRC, alpha, I1, I2, beta1, I3, I4)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GRC(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: alpha, beta1
          Integer, INTENT(IN) :: I1,I2,I3,I4

!           Corr=(dble(N_SUN)*(alpha*Grc(I1,I2,1)+conjg(alpha)*Grc(I1,I2,2))*(beta1*Grc(I3,I4,1)+conjg(beta1)*Grc(I3,I4,2))&
!                       &+(alpha*beta1*Grc(I1,I4,1)*Gr(I2,I3,1)+conjg(alpha*beta1)*Grc(I1,I4,2)*Gr(I2,I3,2)))*dble(N_SUN)
          Corr=(dble(N_SUN)*(alpha*Grc(I1,I2,1))*(beta1*Grc(I3,I4,1))&
                      &+(alpha*beta1*Grc(I1,I4,1)*Gr(I2,I3,1)))*dble(N_SUN)

        end function Corr
!==========================================================   

        Complex (Kind=kind(0.d0)) Function   CorrT(G00, GTT, GT0, G0T, alpha, I1, I2, beta1, I3, I4)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: G00(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: G0T(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: alpha, beta1
          Integer, INTENT(IN) :: I1,I2,I3,I4
          Complex (Kind=kind(0.d0)) :: Z, Delta1, delta2

          Z=dble(N_SUN)
          delta1=0
          if(I1==I2) delta1=1
          delta2=0
          if(I3==I4) delta2=1
!           CorrT=Z*Z*(alpha*(delta1 - GTT(I2,I1,1)) + conjg(alpha)*(delta1 - GTT(I2,I1,2)))*&
!                     &( beta1*(delta2 - G00(I4,I3,1)) + conjg(beta1)*(delta2 - G00(I4,I3,2))) - &
!                     &     Z * ( alpha*beta1*GT0(I2,I3,1)*G0T(I4,I1,1) + conjg(alpha*beta1)*GT0(I2,I3,2)*G0T(I4,I1,2) )
          CorrT=Z*Z*(alpha*(delta1 - GTT(I2,I1,1)) )*&
                    &( beta1*(delta2 - G00(I4,I3,1)) ) - &
                    &     Z * ( alpha*beta1*GT0(I2,I3,1)*G0T(I4,I1,1) )
                    ! Z * ((DeltaI - GTT(I2,I1,1))*(DeltaJ - G00(J2,J1,1)) - GT0(I2,J1,1)*G0T(J2,I1,1)) * ZP* ZS
!           (N_SUN*Grc(I1,I2,1)*Grc(I3,I4,1)+Grc(I1,I4,1)*Gr(I2,I3,1))*N_SUN

        end function CorrT

!==========================================================   

        ! this routine considers ops C^\dag_I1 C_I2 + h.c. and i*C^\dag_I1 C_I2 + h.c. as den and curr
        ! similar I1->I3 and I2->I4 for second Op_set
        ! and returns the respective den-den and curr-curr correlation function
        ! ATTENTION: I'm expoloiting the for nf=2, alpha and beta1 are complex conjugated!!!

        subroutine   Bond_obs(GR, GRC, alpha, I1, I2, beta1, I3, I4, Den, Curr)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN)  :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN)  :: GRC(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: alpha, beta1
          Complex (Kind=Kind(0.d0)), INTENT(OUT) :: Den, Curr
          Integer, INTENT(IN) :: I1,I2,I3,I4
          Complex (Kind=Kind(0.d0)) :: tmp, alpha_tmp, beta1_tmp

          Den=0.d0
          Curr=0.d0
          alpha_tmp=alpha
          beta1_tmp=beta1
          
          tmp=Corr(GR, GRC, alpha_tmp, I1, I2, beta1_tmp, I3, I4)
          Den=Den+tmp
          Curr=Curr-tmp
          
          alpha_tmp=conjg(alpha_tmp)
          tmp=Corr(GR, GRC, alpha_tmp, I2, I1, beta1_tmp, I3, I4)
          Den=Den+tmp
          Curr=Curr+tmp
          
          alpha_tmp=conjg(alpha_tmp)
          beta1_tmp=conjg(beta1_tmp)
          tmp=Corr(GR, GRC, alpha_tmp, I1, I2, beta1_tmp, I4, I3)
          Den=Den+tmp
          Curr=Curr+tmp
          
          alpha_tmp=conjg(alpha_tmp)
          tmp=Corr(GR, GRC, alpha_tmp, I2, I1, beta1_tmp, I4, I3)
          Den=Den+tmp
          Curr=Curr-tmp

        end subroutine Bond_obs
        

!==========================================================   

        ! this routine considers ops C^\dag_I1 C_I2 + h.c. and i*C^\dag_I1 C_I2 + h.c. as den and curr
        ! similar I1->I3 and I2->I4 for second Op_set
        ! and returns the respective den-den and curr-curr correlation function
        ! ATTENTION: I'm expoloiting the for nf=2, alpha and beta1 are complex conjugated!!!

        subroutine   Bond_obsT(G00, GTT, GT0, G0T, alpha, I1, I2, beta1, I3, I4, Den, Curr)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: G00(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: G0T(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: alpha, beta1
          Complex (Kind=Kind(0.d0)), INTENT(OUT) :: Den, Curr
          Integer, INTENT(IN) :: I1,I2,I3,I4
          Complex (Kind=Kind(0.d0)) :: tmp, alpha_tmp, beta1_tmp

          Den=0.d0
          Curr=0.d0
          alpha_tmp=alpha
          beta1_tmp=beta1
          
          tmp=CorrT(G00, GTT, GT0, G0T, alpha_tmp, I1, I2, beta1_tmp, I3, I4)
          Den=Den+tmp
          Curr=Curr-tmp
          
          alpha_tmp=conjg(alpha_tmp)
          tmp=CorrT(G00, GTT, GT0, G0T, alpha_tmp, I2, I1, beta1_tmp, I3, I4)
          Den=Den+tmp
          Curr=Curr+tmp
          
          alpha_tmp=conjg(alpha_tmp)
          beta1_tmp=conjg(beta1_tmp)
          tmp=CorrT(G00, GTT, GT0, G0T, alpha_tmp, I1, I2, beta1_tmp, I4, I3)
          Den=Den+tmp
          Curr=Curr+tmp
          
          alpha_tmp=conjg(alpha_tmp)
          tmp=CorrT(G00, GTT, GT0, G0T, alpha_tmp, I2, I1, beta1_tmp, I4, I3)
          Den=Den+tmp
          Curr=Curr-tmp

        end subroutine Bond_obsT
        


!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)  
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase  
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice 
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, Zn, weight, tmp, Den, Curr, delta
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY, alpha, beta1, flux, Zdo
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n, K, K1, L, I2, J2, I_tmp, nf2
          
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

          ! Compute scalar observables.
          ! last one is meassured in obserT
          Do I = 1,Size(Obs_scal,1)-1
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
             

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do n  = 1,Size(Op_T,1)
             Do nf = 1,N_FL
                Do I = 1,Size(Op_T(n,nf)%O,1)
                   Do J = 1,Size(Op_T(n,nf)%O,2)
                      Zkin = Zkin +  Op_T(n,nf)%O(i, j)*Grc( Op_T(n,nf)%P(I), Op_T(n,nf)%P(J), nf )
                   ENddo
                Enddo
             Enddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS

          Zn=cmplx(dble(N_sun),0.d0,kind(0.d0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do n = 1,size(OP_V,1)-1
                weight=-Op_V(n,nf)%g**2 /dtau
                Do J = 1,Op_V(n,nf)%N
                   J1 = Op_V(n,nf)%P(J)
                   DO I = 1,Op_V(n,nf)%N
                      if (abs(Op_V(n,nf)%O(i,j)) >= 0.00001) then
                      I1 = Op_V(n,nf)%P(I)
                      ZPot  = ZPot  + N_FL*2.d0*Zn*Op_V(n,nf)%alpha*weight*Op_V(n,nf)%O(i,j)*GRC(I1,J1,nf)
                      
                      Do nf2 = 1,N_FL
                        Delta=0.d0
                        if (nf==nf2) Delta=1.d0
                        Do K = 1,Op_V(n,nf)%N
                          K1 = Op_V(n,nf)%P(K)
                          DO L = 1,Op_V(n,nf)%N
                            if (abs(Op_V(n,nf)%O(k,l)) >= 0.00001) then
                            L1 = Op_V(n,nf)%P(L)
                            tmp =  (   delta*GRC(I1,L1,nf) * GR (J1,K1,nf)      +  &
                                  &    Zn * GRC(I1,J1,nf) * GRC(K1,L1,nf2)         )
                            ZPot  = ZPot  + weight*Op_V(n,nf)%O(i,j)*Op_V(n,nf)%O(k,l)*tmp
                            endif
                          Enddo
                        ENddo
                      enddo
                      endif
                   Enddo
                ENddo
!                 if ( .not. Projector) 
                ZPot  = ZPot  + N_FL*weight*(Op_V(n,nf)%alpha**2)*Zn
             Enddo
          Enddo
          Zpot=Zn*Zpot
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS

!           Zdo = cmplx(0.d0,0.d0, kind(0.D0))
!           Do I = 1,Ndim
!             Zdo = Zdo + Grc(i,i,1)*Grc(i,i,2) 
!           enddo
!           Zdo = -2.d0*Zdo* dble(N_SUN)/dble(Ndim)
!           Obs_scal(6)%Obs_vec(1)  =    Obs_scal(6)%Obs_vec(1) + Zrho/dble(Ndim) * ZP*ZS
!           Obs_scal(6)%Obs_vec(2)  =    Obs_scal(6)%Obs_vec(2) + Zdo * ZP*ZS


!           !Kx
!           flux = 1.0d0
!           if (ham_flux) flux = cmplx(1.d0, 1.d0, kind(0.d0))/sqrt(2.d0)
!           Zkin=0.d0
!           Do I=1,Latt%N
!             do no_I=1,N_curr_bonds
!               I2=invlist(I,1)
!               if (no_I>3) I2=invlist(I,2)
!               
!               select case (no_I)
!               case (1)
!                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                   alpha=-ham_t*flux
!               case (2)
!                   I1 = invlist(Latt%nnlist(I,1, 0),1)
!                   alpha=-ham_t2
!               case (3)
!                   I1 = invlist(Latt%nnlist(I,0, 1),1) 
!                   alpha=ham_t2
!               case (4)
!                   I1 = invlist(Latt%nnlist(I,1, 0),1) 
!                   alpha=-ham_t*conjg(flux) 
!               case (5)
!                   I1 = invlist(Latt%nnlist(I,1, 0),2) 
!                   alpha=ham_t2  
!               case (6)
!                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                   alpha=-ham_t2
!               case (7)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,1, 0)
!                   I2=invlist(I,1)
!                   I1 = invlist(I_tmp,1) 
!                   alpha=-ham_t3*4.d0
!               case (8)
!                   I_tmp=Latt%nnlist(I,0, 1)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I2=invlist(I,1)
!                   I1 = invlist(I_tmp,1) 
!                   alpha=-ham_t3*4.d0
!               case (9)
!                   I2=invlist(I,2)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,1, 0)
!                   I1 = invlist(I_tmp,2)
!                   alpha=-ham_t3*4.d0
!               case (10)
!                   I2=invlist(I,2)
!                   I_tmp=Latt%nnlist(I,0, 1)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I1 = invlist(I_tmp,2)
!                   alpha=-ham_t3*4.d0
!               case (11)
!                   I2=invlist(I,1)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I1 = invlist(I_tmp,1)
!                   alpha=-ham_t4*4.d0
!               case (12)
!                   I2=invlist(I,2)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I1 = invlist(I_tmp,2)
!                   alpha=-ham_t4*4.d0
!               case default
!                   Write(6,*) ' Error in  Ham_Hop '  
!                   Stop
!               end select
!               
!               Z = ( alpha*(GRC(I1,I2,1)) + conjg(alpha)*(GRC(I2,I1,1)) ) 
!               Zkin =  Zkin + Z
!             enddo
!           enddo 
!           Zkin = Zkin * dble(N_SUN)/dble(Latt%N)
!           Obs_scal(5)%Obs_vec(1)  =    Obs_scal(5)%Obs_vec(1) + Zkin * ZP*ZS

          
          ! Compute spin-spin, Green, and den-den correlation functions 
          ! ATTENTION last three are sus, controlled in ObserT
          DO I = 1,Size(Obs_eq,1)-3
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          Do I1 = 1,Ndim
            I    = List(I1,1)
            no_I = List(I1,2)
            Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                ! Green
                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                    &                          ( GRC(I1,J1,1) ) *  ZP*ZS 
                !Den
                Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                    & (   GRC(I1,J1,1) * GR(I1,J1,1)    + &
                    &   (GRC(I1,I1,1))*(GRC(J1,J1,1))     ) * ZP*ZS
            enddo
            Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1)) * ZP*ZS
          enddo
               
          
!           flux = 1.0d0
!           if (ham_flux) flux = cmplx(1.d0, 1.d0, kind(0.d0))/sqrt(2.d0)
!           Do I=1,Latt%N
!             Do J=1,Latt%N
!               imj = latt%imj(I,J)
!               do no_I=1,N_curr_bonds
!                 I2=invlist(I,1)
!                 if (no_I>3) I2=invlist(I,2)
!                 
!                 select case (no_I)
!                 case (1)
!                     I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                     alpha=ham_t*flux
!                 case (2)
!                     I1 = invlist(Latt%nnlist(I,1, 0),1)
!                     alpha=ham_t2
!                 case (3)
!                     I1 = invlist(Latt%nnlist(I,0, 1),1) 
!                     alpha=-ham_t2
!                 case (4)
!                     I1 = invlist(Latt%nnlist(I,1, 0),1) 
!                     alpha=ham_t*conjg(flux) 
!                 case (5)
!                     I1 = invlist(Latt%nnlist(I,1, 0),2) 
!                     alpha=-ham_t2  
!                 case (6)
!                     I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                     alpha=+ham_t2
!                 case (7)
!                     I_tmp=Latt%nnlist(I,1, 0)
!                     I_tmp=Latt%nnlist(I_tmp,1, 0)
!                     I2=invlist(I,1)
!                     I1 = invlist(I_tmp,1) 
!                     alpha=ham_t3*2.d0
!                 case (8)
!                     I_tmp=Latt%nnlist(I,0, 1)
!                     I_tmp=Latt%nnlist(I_tmp,0, 1)
!                     I2=invlist(I,1)
!                     I1 = invlist(I_tmp,1) 
!                     alpha=ham_t3*2.d0
!                 case (9)
!                     I2=invlist(I,2)
!                     I_tmp=Latt%nnlist(I,1, 0)
!                     I_tmp=Latt%nnlist(I_tmp,1, 0)
!                     I1 = invlist(I_tmp,2)
!                     alpha=ham_t3*2.d0
!                 case (10)
!                     I2=invlist(I,2)
!                     I_tmp=Latt%nnlist(I,0, 1)
!                     I_tmp=Latt%nnlist(I_tmp,0, 1)
!                     I1 = invlist(I_tmp,2)
!                     alpha=ham_t3*2.d0
!                 case (11)
!                     I2=invlist(I,1)
!                     I_tmp=Latt%nnlist(I,1, 0)
!                     I_tmp=Latt%nnlist(I_tmp,0, 1)
!                     I1 = invlist(I_tmp,1)
!                     alpha=ham_t4*2.d0
!                 case (12)
!                     I2=invlist(I,2)
!                     I_tmp=Latt%nnlist(I,1, 0)
!                     I_tmp=Latt%nnlist(I_tmp,0, 1)
!                     I1 = invlist(I_tmp,2)
!                     alpha=ham_t4*2.d0
!                 case default
!                     Write(6,*) ' Error in  Ham_Hop '  
!                     Stop
!                 end select
!                 
!                 do no_j=1,N_curr_bonds
!                   J2=invlist(J,1)
!                   if (no_J>3) J2=invlist(J,2)
!                 
!                   select case (no_J)
!                   case (1)
!                       J1 = invlist(Latt%nnlist(J,0, 1),2)
!                       beta1=ham_t*flux 
!                   case (2)
!                       J1 = invlist(Latt%nnlist(J,1, 0),1)
!                       beta1=ham_t2
!                   case (3)
!                       J1 = invlist(Latt%nnlist(J,0, 1),1) 
!                       beta1=-ham_t2
!                   case (4)
!                       J1 = invlist(Latt%nnlist(J,1, 0),1) 
!                       beta1=ham_t*conjg(flux)
!                   case (5)
!                       J1 = invlist(Latt%nnlist(J,1, 0),2)  
!                       beta1=-ham_t2
!                   case (6)
!                       J1 = invlist(Latt%nnlist(J,0, 1),2)
!                       beta1=ham_t2
!                   case (7)
!                       I_tmp=Latt%nnlist(J,1, 0)
!                       I_tmp=Latt%nnlist(I_tmp,1, 0)
!                       J2=invlist(J,1)
!                       J1 = invlist(I_tmp,1) 
!                       beta1=ham_t3*2.d0
!                   case (8)
!                       I_tmp=Latt%nnlist(J,0, 1)
!                       I_tmp=Latt%nnlist(I_tmp,0, 1)
!                       J2=invlist(J,1)
!                       J1 = invlist(I_tmp,1) 
!                       beta1=ham_t3*2.d0
!                   case (9)
!                       J2=invlist(J,2)
!                       I_tmp=Latt%nnlist(J,1, 0)
!                       I_tmp=Latt%nnlist(I_tmp,1, 0)
!                       J1 = invlist(I_tmp,2)
!                       beta1=ham_t3*2.d0
!                   case (10)
!                       J2=invlist(J,2)
!                       I_tmp=Latt%nnlist(J,0, 1)
!                       I_tmp=Latt%nnlist(I_tmp,0, 1)
!                       J1 = invlist(I_tmp,2)
!                       beta1=ham_t3*2.d0
!                   case (11)
!                       J2=invlist(J,1)
!                       I_tmp=Latt%nnlist(J,1, 0)
!                       I_tmp=Latt%nnlist(I_tmp,0, 1)
!                       J1 = invlist(I_tmp,1)
!                       beta1=ham_t4*2.d0
!                   case (12)
!                       J2=invlist(J,2)
!                       I_tmp=Latt%nnlist(J,1, 0)
!                       I_tmp=Latt%nnlist(I_tmp,0, 1)
!                       J1 = invlist(I_tmp,2)
!                       beta1=ham_t4*2.d0
!                   case default
!                       Write(6,*) ' Error in  Ham_Hop '  
!                       Stop
!                   end select
!                   
!                   call Bond_obs(GR, GRC, alpha, I1, I2, beta1, J1, J2, Den, Curr)
!                   Z=Curr
!                   Obs_eq(7)%Obs_Latt(imj,1,no_i,no_j) =  Obs_eq(7)%Obs_Latt(imj,1,no_i,no_j) + Z * ZP * ZS
!                 enddo
!               enddo
!             enddo
! !             do no_I=1,N_curr_bonds
! !               I2=invlist(I,1)
! !               if (no_I>3) I2=invlist(I,2)
! !               
! !               select case (no_I)
! !               case (1)
! !                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
! !                   alpha=ham_t*flux
! !               case (2)
! !                   I1 = invlist(Latt%nnlist(I,1, 0),1)
! !                   alpha=ham_t2
! !               case (3)
! !                   I1 = invlist(Latt%nnlist(I,0, 1),1) 
! !                   alpha=-ham_t2
! !               case (4)
! !                   I1 = invlist(Latt%nnlist(I,1, 0),1) 
! !                   alpha=ham_t*conjg(flux) 
! !               case (5)
! !                   I1 = invlist(Latt%nnlist(I,1, 0),2) 
! !                   alpha=-ham_t2  
! !               case (6)
! !                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
! !                   alpha=+ham_t2
! !               case (7)
! !                   I_tmp=Latt%nnlist(I,1, 0)
! !                   I_tmp=Latt%nnlist(I_tmp,1, 0)
! !                   I2=invlist(I,1)
! !                   I1 = invlist(I_tmp,1) 
! !                   alpha=-ham_t3
! !               case (8)
! !                   I_tmp=Latt%nnlist(I,0, 1)
! !                   I_tmp=Latt%nnlist(I_tmp,0, 1)
! !                   I2=invlist(I,1)
! !                   I1 = invlist(I_tmp,1) 
! !                   alpha=-ham_t3
! !               case (9)
! !                   I2=invlist(I,2)
! !                   I_tmp=Latt%nnlist(I,1, 0)
! !                   I_tmp=Latt%nnlist(I_tmp,1, 0)
! !                   I1 = invlist(I_tmp,2)
! !                   alpha=-ham_t3
! !               case (10)
! !                   I2=invlist(I,2)
! !                   I_tmp=Latt%nnlist(I,0, 1)
! !                   I_tmp=Latt%nnlist(I_tmp,0, 1)
! !                   I1 = invlist(I_tmp,2)
! !                   alpha=-ham_t3
! !               case default
! !                   Write(6,*) ' Error in  Ham_Hop '  
! !                   Stop
! !               end select
! !               
! !               Z = cmplx(0.d0,1.d0,kind(0.d0))*&
! !                   &(alpha*(GRC(I1,I2,1) - GRC(I2,I1,2)) + conjg(alpha)*(GRC(I1,I2,2) - GRC(I2,I1,1))) 
! !               Obs_eq(7)%Obs_Latt0(no_I) =  Obs_eq(7)%Obs_Latt0(no_I) + Z * ZP * ZS
! !             enddo
!           enddo 
!           
!           alpha=1.d0
!           beta1=1.d0
!           Do I=1,Latt%N
!             Do J=1,Latt%N
!               imj = latt%imj(I,J)
!               I2 = Invlist(I,1)
!               I1 = I2
!               Do no_I = 1,Latt_unit%N_coord
!                 select case (no_I)
!                 case (1)
!                     I1 = invlist(I,2) 
!                 case (2)
!                     I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                 case (3)
!                     I1 = invlist(Latt%nnlist(I,-1,1),2) 
!                 case (4)
!                     I1 = invlist(Latt%nnlist(I,-1,0),2) 
!                 case default
!                     Write(6,*) ' Error in  Ham_Hop '  
!                     Stop
!                 end select
!                 
!                 J2=invlist(J,1)
!                 J1 = J2
!                 do no_j=1,Latt_unit%N_coord
!                   select case (no_J)
!                   case (1)
!                       J1 = invlist(J,2) 
!                   case (2)
!                       J1 = invlist(Latt%nnlist(J,0, 1),2) 
!                   case (3)
!                       J1 = invlist(Latt%nnlist(J,-1,1),2) 
!                   case (4)
!                       J1 = invlist(Latt%nnlist(J,-1,0),2) 
!                   case default
!                       Write(6,*) ' Error in  Ham_Hop '  
!                       Stop
!                   end select
!                   
!                   call Bond_obs(GR, GRC, alpha, I1, I2, beta1, J1, J2, Den, Curr)
!                   Z=Den
!                   Obs_eq(8)%Obs_Latt(imj,1,no_i,no_j) =  Obs_eq(8)%Obs_Latt(imj,1,no_i,no_j) + Z * ZP * ZS
!                 enddo
!               enddo
!             enddo
!             
!             I2 = Invlist(I,1)
!             I1 = I2
!             do no_I=1,Latt_unit%N_coord
!               select case (no_I)
!               case (1)
!                   I1 = invlist(I,2) 
!               case (2)
!                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
!               case (3)
!                   I1 = invlist(Latt%nnlist(I,-1,1),2) 
!               case (4)
!                   I1 = invlist(Latt%nnlist(I,-1,0),2) 
!               case default
!                   Write(6,*) ' Error in  Ham_Hop '  
!                   Stop
!               end select
!               
!               Z = (alpha*(GRC(I1,I2,1) + GRC(I2,I1,2)) + conjg(alpha)*(GRC(I1,I2,2) + GRC(I2,I1,1))) 
!               Obs_eq(8)%Obs_Latt0(no_I) =  Obs_eq(8)%Obs_Latt0(no_I) + Z * ZP * ZS
!             enddo
!           enddo 

        end Subroutine Obser
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)  
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )> 
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)> 
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )> 
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)> 
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase  
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, Den, Curr, alpha, beta1, flux, weight, drude
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, I2, J2, sus_ind, I_tmp
          real(Kind=kind(0.d0)), parameter :: pi=acos(-1.d0)
          real(Kind=kind(0.d0)) :: ang

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          sus_ind = size(Obs_eq,1)
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
             DO I = Size(Obs_eq,1)-2,Size(Obs_eq,1)
                Obs_eq(I)%N        = Obs_eq(I)%N + 1
                Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
             ENDDO
             Obs_scal(size(Obs_scal,1))%N         =  Obs_scal(size(Obs_scal,1))%N + 1
             Obs_scal(size(Obs_scal,1))%Ave_sign  =  Obs_scal(size(Obs_scal,1))%Ave_sign + Real(ZS,kind(0.d0))
          endif
          
          weight=dtau
          if ((NT == 0 .or. NT == Ltrot) .and. .not. projector) weight=0.5d0*dtau
          Do I1 = 1,Ndim
            I    = List(I1,1)
            no_I = List(I1,2)
            Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                !Green
                Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS
                !Green_up
                Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) + GT0(I1,J1,1) * ZP* ZS
                !Green_down
                Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) + GT0(I1,J1,2) * ZP* ZS
                !Green_0T (only for projector)
                if (projector) then
                  Obs_tau(size(Obs_tau,1))%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(size(Obs_tau,1))%Obs_Latt(imj,nt+1,no_I,no_J)  &
                      &   -   ( G0T(I1,J1,1) + G0T(I1,J1,2) ) * ZP* ZS
                endif

                !SpinZ
                Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  ( &
                    &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                    &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS
!                 !SpinXY
!                 Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                     &  - &
!                     &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                !Den
                Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  (                                        &  
                    &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                    &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                    & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS    
                !Den_sus
                Obs_eq(sus_ind-2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(sus_ind-2)%Obs_Latt(imj,1,no_I,no_J)  &
                    & +  weight*(                                        &  
                    &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                    &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                    & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS   
                !Spin_sus
                Obs_eq(sus_ind-1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(sus_ind-1)%Obs_Latt(imj,1,no_I,no_J)  &
                    & -  weight*( G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1) )*ZP*ZS
                !SC
                Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & +  G0T(J1,I1,1) * G0T(J1,I1,2) *ZP*ZS    
                !SC_sym
                Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                    & + 0.25d0*( G0T(J1,I1,1) * G0T(J1,I1,2) + GT0(I1,J1,1) * GT0(I1,J1,2) )*ZP*ZS   

            enddo
          
            Obs_tau(7)%Obs_Latt0(no_I) =  Obs_tau(7)%Obs_Latt0(no_I) + &
                  &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
            Obs_eq(sus_ind-2)%Obs_Latt0(no_I) =  Obs_eq(sus_ind-2)%Obs_Latt0(no_I) + &
                  &       weight/sqrt(beta)*(cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
          Enddo
          
          weight=dtau
          drude=0.d0
          if ((NT == 0 .or. NT == Ltrot) .and. .not. projector) weight=0.5d0*dtau
          flux = 1.0d0
          if (ham_flux) flux = cmplx(1.d0, 1.d0, kind(0.d0))/sqrt(2.d0)
          Do I=1,Latt%N
            Do J=1,Latt%N
              imj = latt%imj(I,J)
              do no_I=1,N_curr_bonds
                I2=invlist(I,1)
                if (no_I>3) I2=invlist(I,2)
                
                select case (no_I)
                case (1)
                    I1 = invlist(Latt%nnlist(I,0, 1),2) 
                    alpha=ham_t*flux
                case (2)
                    I1 = invlist(Latt%nnlist(I,1, 0),1)
                    alpha=ham_t2
                case (3)
                    I1 = invlist(Latt%nnlist(I,0, 1),1) 
                    alpha=-ham_t2
                case (4)
                    I1 = invlist(Latt%nnlist(I,1, 0),1) 
                    alpha=ham_t*conjg(flux) 
                case (5)
                    I1 = invlist(Latt%nnlist(I,1, 0),2) 
                    alpha=-ham_t2  
                case (6)
                    I1 = invlist(Latt%nnlist(I,0, 1),2) 
                    alpha=+ham_t2
                case (7)
                    I_tmp=Latt%nnlist(I,1, 0)
                    I_tmp=Latt%nnlist(I_tmp,1, 0)
                    I2=invlist(I,1)
                    I1 = invlist(I_tmp,1) 
                    alpha=ham_t3*2.d0
                case (8)
                    I_tmp=Latt%nnlist(I,0, 1)
                    I_tmp=Latt%nnlist(I_tmp,0, 1)
                    I2=invlist(I,1)
                    I1 = invlist(I_tmp,1) 
                    alpha=ham_t3*2.d0
                case (9)
                    I2=invlist(I,2)
                    I_tmp=Latt%nnlist(I,1, 0)
                    I_tmp=Latt%nnlist(I_tmp,1, 0)
                    I1 = invlist(I_tmp,2)
                    alpha=ham_t3*2.d0
                case (10)
                    I2=invlist(I,2)
                    I_tmp=Latt%nnlist(I,0, 1)
                    I_tmp=Latt%nnlist(I_tmp,0, 1)
                    I1 = invlist(I_tmp,2)
                    alpha=ham_t3*2.d0
                case (11)
                    I2=invlist(I,1)
                    I_tmp=Latt%nnlist(I,1, 0)
                    I_tmp=Latt%nnlist(I_tmp,0, 1)
                    I1 = invlist(I_tmp,1)
                    alpha=ham_t4*2.d0
                case (12)
                    I2=invlist(I,2)
                    I_tmp=Latt%nnlist(I,1, 0)
                    I_tmp=Latt%nnlist(I_tmp,0, 1)
                    I1 = invlist(I_tmp,2)
                    alpha=ham_t4*2.d0
                case default
                    Write(6,*) ' Error in  Ham_Hop '  
                    Stop
                end select
                
                do no_j=1,N_curr_bonds
                  J2=invlist(J,1)
                  if (no_J>3) J2=invlist(J,2)
                
                  select case (no_J)
                  case (1)
                      J1 = invlist(Latt%nnlist(J,0, 1),2)
                      beta1=ham_t*flux 
                  case (2)
                      J1 = invlist(Latt%nnlist(J,1, 0),1)
                      beta1=ham_t2
                  case (3)
                      J1 = invlist(Latt%nnlist(J,0, 1),1) 
                      beta1=-ham_t2
                  case (4)
                      J1 = invlist(Latt%nnlist(J,1, 0),1) 
                      beta1=ham_t*conjg(flux)
                  case (5)
                      J1 = invlist(Latt%nnlist(J,1, 0),2)  
                      beta1=-ham_t2
                  case (6)
                      J1 = invlist(Latt%nnlist(J,0, 1),2)
                      beta1=ham_t2
                  case (7)
                      I_tmp=Latt%nnlist(J,1, 0)
                      I_tmp=Latt%nnlist(I_tmp,1, 0)
                      J2=invlist(J,1)
                      J1 = invlist(I_tmp,1) 
                      beta1=ham_t3*2.d0
                  case (8)
                      I_tmp=Latt%nnlist(J,0, 1)
                      I_tmp=Latt%nnlist(I_tmp,0, 1)
                      J2=invlist(J,1)
                      J1 = invlist(I_tmp,1) 
                      beta1=ham_t3*2.d0
                  case (9)
                      J2=invlist(J,2)
                      I_tmp=Latt%nnlist(J,1, 0)
                      I_tmp=Latt%nnlist(I_tmp,1, 0)
                      J1 = invlist(I_tmp,2)
                      beta1=ham_t3*2.d0
                  case (10)
                      J2=invlist(J,2)
                      I_tmp=Latt%nnlist(J,0, 1)
                      I_tmp=Latt%nnlist(I_tmp,0, 1)
                      J1 = invlist(I_tmp,2)
                      beta1=ham_t3*2.d0
                  case (11)
                      J2=invlist(J,1)
                      I_tmp=Latt%nnlist(J,1, 0)
                      I_tmp=Latt%nnlist(I_tmp,0, 1)
                      J1 = invlist(I_tmp,1)
                      beta1=ham_t4*2.d0
                  case (12)
                      J2=invlist(J,2)
                      I_tmp=Latt%nnlist(J,1, 0)
                      I_tmp=Latt%nnlist(I_tmp,0, 1)
                      J1 = invlist(I_tmp,2)
                      beta1=ham_t4*2.d0
                  case default
                      Write(6,*) ' Error in  Ham_Hop '  
                      Stop
                  end select
                  
                  
                  call Bond_obsT(G00, GTT, GT0, G0T, alpha, I1, I2, beta1, J1, J2, Den, Curr)
                  Z=Curr
                  drude=drude + Curr
                  Obs_tau(3)%Obs_Latt(imj,nt+1,no_i,no_j) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_i,no_j) + Z * ZP * ZS
                  Obs_eq(sus_ind)%Obs_Latt(imj,1,no_i,no_j) =  Obs_eq(sus_ind)%Obs_Latt(imj,1,no_i,no_j) + weight * Z * ZP * ZS
                enddo
              enddo
            enddo
!             do no_I=1,N_curr_bonds
!               I2=invlist(I,1)
!               if (no_I>3) I2=invlist(I,2)
!               
!               select case (no_I)
!               case (1)
!                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                   alpha=ham_t*flux
!               case (2)
!                   I1 = invlist(Latt%nnlist(I,1, 0),1)
!                   alpha=ham_t2
!               case (3)
!                   I1 = invlist(Latt%nnlist(I,0, 1),1) 
!                   alpha=-ham_t2
!               case (4)
!                   I1 = invlist(Latt%nnlist(I,1, 0),1) 
!                   alpha=ham_t*conjg(flux) 
!               case (5)
!                   I1 = invlist(Latt%nnlist(I,1, 0),2) 
!                   alpha=-ham_t2  
!               case (6)
!                   I1 = invlist(Latt%nnlist(I,0, 1),2) 
!                   alpha=+ham_t2
!               case (7)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,1, 0)
!                   I2=invlist(I,1)
!                   I1 = invlist(I_tmp,1) 
!                   alpha=-ham_t3
!               case (8)
!                   I_tmp=Latt%nnlist(I,0, 1)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I2=invlist(I,1)
!                   I1 = invlist(I_tmp,1) 
!                   alpha=-ham_t3
!               case (9)
!                   I2=invlist(I,2)
!                   I_tmp=Latt%nnlist(I,1, 0)
!                   I_tmp=Latt%nnlist(I_tmp,1, 0)
!                   I1 = invlist(I_tmp,2)
!                   alpha=-ham_t3
!               case (10)
!                   I2=invlist(I,2)
!                   I_tmp=Latt%nnlist(I,0, 1)
!                   I_tmp=Latt%nnlist(I_tmp,0, 1)
!                   I1 = invlist(I_tmp,2)
!                   alpha=-ham_t3
!               case default
!                   Write(6,*) ' Error in  Ham_Hop '  
!                   Stop
!               end select
!               
!               Z = - cmplx(0.d0,1.d0,kind(0.d0))*&
!                   &(alpha*(GTT(I1,I2,1) - GTT(I2,I1,2)) + conjg(alpha)*(GTT(I1,I2,2) - GTT(I2,I1,1)) ) 
!               Obs_tau(3)%Obs_Latt0(no_I) =  Obs_tau(3)%Obs_Latt0(no_I) + Z * ZP * ZS
!               Obs_eq(sus_ind)%Obs_Latt0(no_I) =  Obs_eq(sus_ind)%Obs_Latt0(no_I) + weight* Z * ZP * ZS
!             enddo
          enddo
          do I=0,size(Obs_scal(size(Obs_scal,1))%Obs_vec)-1
            ang=nt*2.d0*pi/dble(Ltrot)
            if (projector) ang=nt*pi/dble(Ltrot-2*Thtrot)
            Obs_scal(size(Obs_scal,1))%Obs_vec(I+1)=Obs_scal(size(Obs_scal,1))%Obs_vec(I+1)&
                                    &+cmplx(cos(I*ang),sin(I*ang),kind(0.d0))*weight*drude/dble(Latt%N)*ZP*ZS
          enddo
          
        end Subroutine OBSERT
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> Prints out the bins.  No need to change this routine.
!-------------------------------------------------------------------
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau
          
          !Local 
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I),Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I),Latt,dtau,Group_Comm)
          enddo
          If (Ltau  == 1 ) then
            Do I = 1,Size(Obs_tau,1)
              if(I==3) then !compactify for curr to save disc space
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau,Group_Comm,Bond_pos)
              else
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau,Group_Comm)
              endif
            enddo
          endif

        end Subroutine Pr_obs

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> Initializes observables to zero before each bins.  No need to change
!> this routine.
!-------------------------------------------------------------------
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          ! Local 
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Specify a global move on a given time slice tau.
!>
!> @details
!> @param[in] ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!> @param[out] T0_Proposal_ratio, Real
!> \verbatim
!>  T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)
!> \endverbatim
!> @param[out] S0_ratio, Real
!> \verbatim
!>  S0_ratio = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> \endverbatim
!> @param[out] Flip_length  Integer
!> \verbatim
!>  Number of flips stored in the first  Flip_length entries of the array Flip_values.
!>  Has to be smaller than NDIM
!> \endverbatim
!> @param[out] Flip_list  Integer(Ndim)
!> \verbatim
!>  List of spins to be flipped: nsigma%f(Flip_list(1),ntau) ... nsigma%f(Flip_list(Flip_Length),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!> @param[out] Flip_value  Real(Ndim)
!> \verbatim
!>  Flip_value(:)= nsigma%flip(Flip_list(:),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!--------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)

          
          Implicit none 
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio,  S0_ratio
          Integer                , INTENT(OUT) :: Flip_list(:)
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal

          If (Model == "LRC" ) then
             Call LRC_draw_field(Percent_change, Dtau, nsigma%f(:,ntau), Flip_value,N_SUN)
             Do n = 1,Ndim
                Flip_list(n) = n
                !Write(6,*) Flip_value(n), nsigma%f(n,ntau)
             Enddo
             !Write(6,*)
             Flip_length    = Ndim
             ! T0_Proposal_ration exactly cancels S0_ratio
             T0_Proposal_ratio = 1.d0
             S0_ratio          = 1.d0
          else
             Flip_length = nranf(4)
             do n = 1,flip_length 
                n_op = nranf(size(OP_V,1))
                Flip_list(n)  = n_op
                Flip_value(n) = nsigma%flip(n_op,ntau)
                If ( OP_V(n_op,1)%type == 1 ) then 
                   S0_ratio          =   S0(n_op,ntau,Flip_value(n))
                   T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0_ratio) ! No move prob
                   If ( T0_Proposal > Ranf_wrap() ) then
                      T0_Proposal_ratio =  1.d0 / S0_ratio
                   else
                      T0_Proposal_ratio = 0.d0
                   endif
                else
                   T0_Proposal_ratio = 1.d0
                   S0_ratio          = 1.d0
                endif
             Enddo
          endif
          
        end Subroutine Global_move_tau

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> The user can set the initial field.
!>
!> @details
!> @param[OUT] Initial_field Real(:,:)
!> \verbatim
!>  Upon entry Initial_field is not allocated. If alloacted then it will contain the
!>  the initial field
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine  Hamiltonian_set_nsigma(Initial_field) 
        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field

        
      end Subroutine Hamiltonian_set_nsigma
!--------------------------------------------------------------------

        
      end Module Hamiltonian
