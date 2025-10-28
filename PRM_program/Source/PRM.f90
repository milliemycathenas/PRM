module ProcessVariables

    use Constants
    
    implicit none
  
    integer :: count, index, label, mark
    integer :: row, column
    integer :: BeginPoint, EndPoint
    integer :: bra, ket !count for matrix elemets, also represents bra and ket vector
    integer :: PhiBra, PhiKet
    integer :: CountMu, CountNu
    integer :: countr_max !dimension of rotor Hamiltonian (R = I-J)
    integer :: serialnum !correspond the state and particle number, in place of sequence
    integer, allocatable, dimension(:,:) :: NPassArray, IPassArray
    !storing if a state is catagorized in normal or invesre state for each j shell
    integer :: IfQuasi !if using BCS transformation or not (yes=1, no=0)
    integer :: IfEMTrans !if calculating electromagnetic transitions or not (yes=1, no=0)
    integer :: IfTrunc !if truancation

    real :: KBra, KKet
    real :: SignPass, PassMark !assignment passing
    real :: PSignPass, NSignPass !same as SignPass for proton and neutron respectively

    character(len=50) :: VecFileName !storage of single particle arrangement
    character(len=50) :: SCFileName !storage of strong coupling vector serial number
    character(len=50) :: SPEFileName !judgement for exporting to proton or neutron energy file
    character(len=50) :: SPSFileName !judgement for exporting to proton or neutron state file
    character(len=50) :: RFileName !storage of rotation hamiltion for each I
    character(len=50) :: TFileName !storage of total hamiltion for each I

    integer, allocatable, dimension(:) :: PorN !marking each orbital(tau3) is for proton(1) or neutron(0)
    integer(kind=8), allocatable, dimension(:) :: countsp_max!dimension of single particle Hamiltonian

    real, allocatable, dimension(:,:,:) :: EAlign !reordered single particle energies
    real, allocatable, dimension(:,:,:,:) :: CAlign !reordered single particle states
    
end module ProcessVariables

module Constants

    implicit none

    !program parameters
    integer, parameter :: DIM=1000
    integer, parameter :: DIMA=100000 !dimension for single particle sequence array 
    integer, parameter :: SPoutput=14 !unit of output file of single particle energy
    integer, parameter :: Routput=15 !unit of output file of rotor energy
    integer, parameter :: PSoutput=16 !unit of output file of particle occupation
    integer, parameter :: Ioutput=17 !unit of output file of intrinsic part
    integer, parameter :: Voutput=18 !unit of output file of valence part
    integer, parameter :: Toutput=19 !unit of output file of total part

   
    !math parameters
    real, parameter :: pi=3.1415926

end module Constants

module PhysicalQuantities

    use Constants
    use ProcessVariables

    implicit none

    !!!input parameters begin!!!
        integer :: Mass !mass number of the calculated nucleus
        real :: Omega !magnetic quantum number
        real :: hbaromega_0 !harmonic oscillator energy unit
        real :: NucRadius !radius of the calculated nucleus
        real :: NucProtonNum !proton number(charge number) of the calculated nucleus
        real :: CoriolisF !Coriolis attenuation Factor
        character(len=8) :: Nuclide
        character(len=50) :: NuclideMass
        character(len=50) :: IName
    !!!input parameters end!!!

    !!! single particle energy variables begin !!!
        integer :: J1State, J2State, J3State, J4State !loop index for particle state organization
        integer :: tau3_max !total orbital number of jsp1, jsp2, jsp3...
        integer :: tcount_max !the orbital(tau3) that possesses largest countsp space
        integer(kind=8) :: SPDimNum !single particle matrix dimension number only for normal states

        real :: MOI_0 !moment of inertia J_0
        real :: jsp3_squared !z-axis projection square of total angular momentum
        real :: jsp_squared !square of total angular momentum
        real :: jspup_squared, jspdown_squared !square of angular shift operators
        real :: gamma, gamma_rad !quadrupole deformation parameter gamma, rad is its radian value
        real :: beta, varepsilon !total deformation parameter beta & varepsilon 
        
        integer, allocatable, dimension(:) :: ParOrHole !particle on orbital>ParOrHole(tau3)=1, hole on orbital>ParOrHole(tau3)=-1
        integer, allocatable, dimension(:) :: N !principal quantum number

        real, allocatable, dimension(:,:,:) :: NilVector 
        !Nilsson basis as N-l-j-Omega; NilVector(1,DIM):j value; NilVector(2,DIM):Omega value
        real(kind=4), allocatable, dimension(:,:) :: SPHamiltonian !single particle Hamiltonian matrix element
        real, allocatable, dimension(:) :: DeformC !deformation parameter C=\varepsilon_2(N+3/2)\hbar\omega_0\frac{1}{j(j+1)}
        real, allocatable, dimension(:) :: jsp !total angular momentum quantum number
        real(kind=4), allocatable, dimension(:) :: SPEnergy !energy of the single particle levels
        real, allocatable, dimension(:,:) :: SPState !eigen states pf the single particle levels
    !!! single particle energy variables end !!!

    !!! rotor energy variables begin !!!
        real :: I, IStart, IEnd !total angular momentum of nucleus
        real :: K !magnetic quantum number 
        real :: Irotor3 !projection of I
        real :: Irotor_squared !squared projection of I
        real :: Irotorup_squared, Irotordown_squared ! shift operators of I

        real, allocatable, dimension(:) :: MOI !moment of inertia on x, y, z directions
        real, allocatable, dimension(:) :: A !A1, A2, A3 in H_coll 
        real(kind=4), allocatable, dimension(:,:) :: RHamiltonian !rotor Hamiltonian matrix element
        real(kind=4), allocatable, dimension(:) :: REnergy
        real, allocatable, dimension(:,:) :: RState
    !!! rotor energy variables end !!!

    !!! intrinsic energy variables begin !!!
        integer :: BinaryNumber !number of possible states for the whole valence space
        integer :: TotalStates !number of particle states for normal phi
        integer :: sequence !sequence number of nucleon, proton in front, neutron behind
        integer :: D2Symmetry !value of K-0.5(z1-z2)-0.5(n1-n2)
        integer(kind=8) :: DimensionNumber !dimension of strong coupling basis |IMK$\phi$>

        real :: E_cut !energy for truncation
        real :: E_MPC !MPC energy in multi-particle situation
        real :: E_min !lowest MPC energy in multi-particle situation
        real :: E_difference !E_MPC - E_min

        integer, allocatable, dimension(:) :: StateRemained !state remained after considering proton and neutron number
        integer, allocatable, dimension(:) :: StateSequence !serial number of sorted single particle states
        integer, allocatable, dimension(:) :: StateNormal !number of selected normal states for R2 symmetry
        integer, allocatable, dimension(:) :: Combination!result of C_{2jp+1} and C_{2jn+1}
        integer, allocatable, dimension(:) :: Tau3PNum !valence particle number on each orbit
        integer, allocatable, dimension(:) :: DimensionArray !storing different dimension number for different I
        integer, allocatable, dimension(:) :: KNumber !number of K values conform to R3 symmetry
        integer, allocatable, dimension(:,:) :: JShellIndex
        !combination of different particle states for each single j shell
        integer, allocatable, dimension(:,:) :: NStateCoupled, IStateCoupled 
        !state serial number remained after considering R2 symmetry for each j shell

        integer, allocatable, dimension(:,:,:) :: SPVector !single particle states divided into different j shells for normal states
        integer, allocatable, dimension(:,:,:) :: ISPVector !single particle states divided into different j shells for inverse states
        integer, allocatable, dimension(:,:) :: PhiVector !normal particle states for all single j shells in one vector
        integer, allocatable, dimension(:,:) :: IPhiVector !inverse particle states for all single j shells in one vector
        integer, allocatable, dimension(:,:) :: PNumN, PNumI !proton number on normal or inverse state nu
        
        integer, allocatable, dimension(:) :: ParticleVector !sorted particle states for both neutron and proton&
        !colunm represents for specific state, row represents for sequence
        
        real, allocatable, dimension(:) :: IEnergy !energy of the intrinsic part
        real, allocatable, dimension(:,:) :: IState !eigen state of the intrinsic part
        real(kind=4), allocatable, dimension(:,:) :: IHamiltonian !Intrinsic Hamiltonian matrix elemets
        real, allocatable, dimension(:,:) :: NSCVector !contains allowable K value with corresponding state sequence &
        !NSCVector(1,bra)=K, NSCVector(2,bra)=sequence
        real, allocatable, dimension(:,:) :: ISCVector !contains allowable K value with corresponding inverse state sequence &
        !ISCVector(1,bra)=K, SCVector(2,bra)=sequence that makes SCVector(2,bra) and ISCVector(2,bra) satisfy \delta{\phi\phi'}
    !!! intrinsic energy variables end !!!

    !!! decoupling factor variables begin !!!
        integer :: MuSequence, NuSequence
        !serial number for \chi_mu and \chi_nu, value from 1 to 2*jsp+1 &
        !MuSeApo and NuSeApo stands for \chi_mu' and \chi_nu'
        integer :: MuSeApo, NuSeApo, MuSerial, NuSerial

        real, allocatable, dimension(:,:) :: DecVector !j and Omega in Decoupling factors, start from -jsp-0.5 &
        !DecVector(1,count)=j, DecVector(2,count)=Omega
        real, allocatable, dimension(:,:) :: DecoupleAB, DecoupleC, DecoupleD, DecoupleE, DecoupleF, DecoupleG 
        real, allocatable, dimension(:,:,:) :: C_nu !reordered calculated single particle expansion coefficient of state mu &
        !row=\mu or \nu; column=\Omega
        real, allocatable, dimension(:,:,:) :: DecoupleABSum, DecoupleCSum, DecoupleDSum, DecoupleESum, DecoupleFSum, DecoupleGSum 
        !final result of decouling factors for proton with only positive array index
    !!! decoupling factor variables end !!!

    !!!valence energy variables begin!!!
        integer :: tau3, tau3Apo, orbit !serial number of different orbits
        integer :: MuTau3, NuTau3, MuTau3Apo, NuTau3Apo
        !orbital that the return value of Mu and Nu
        integer :: AnnihMu, AnnihNu, AnnihMuApo, AnnihNuApo

        real :: CoeffNtoI !the coefficient of (-1)**(I-K)(-1)**(...) when change normal state to inverse state
        real :: CoeffAnnih !factor "-1" after operating annihilate operator on single particle state

        integer, allocatable, dimension(:,:) :: BSPVector, KSPVector
        !temporary array storing single particle vectors operated by annihilate operators &
        !in monomer operator case

        real :: NJSMono, NJSPMono, NJSTB
        !hamiltonian elements for Jsquared; JSPMono stands for J squared product of monomer operator
        real :: NJ3SMono, NJ3SPMono, NJ3STB
        !hamiltonian elements for J3squared; J3SPMono stands for J3 squared product of monomer operator
        real :: NJ3Mono
        !hamiltonian elements for J3; J3PMono stands for J3 product of monomer operator
        real :: NJuJdMono, NJuJdPMono, NJuJdTB
        !hamiltonian elements for JupJdown; JuJdPMono stands for JupJdown product of monomer operator
        real :: NJdJuMono, NJdJuPMono, NJdJuTB
        !hamiltonian elemets for JdownJup; JuJdPMono stands for JdownJup product of monomer operator
        real :: NJuMono
        !hamiltonian elemets for Jup; JuTB stands for Jup two-body matrix operator
        real :: NJdMono
        !hamiltonian elemets for Jdown; JdTB stands for Jdown two-body matrix operator

        real :: IJSMono, IJSPMono, IJSTB
        !hamiltonian elements for Jsquared; JSPMono stands for J squared product of monomer operator
        real :: IJ3SMono, IJ3SPMono, IJ3STB
        !hamiltonian elements for J3squared; J3SPMono stands for J3 squared product of monomer operator
        real :: IJ3Mono
        !hamiltonian elements for J3; J3PMono stands for J3 product of monomer operator
        real :: IJuJdMono, IJuJdPMono, IJuJdTB
        !hamiltonian elements for JupJdown; JuJdPMono stands for JupJdown product of monomer operator
        real :: IJdJuMono, IJdJuPMono, IJdJuTB
        !hamiltonian elemets for JdownJup; JuJdPMono stands for JdownJup product of monomer operator
        real :: IJuMono
        !hamiltonian elemets for Jup; JuTB stands for Jup two-body matrix operator
        real :: IJdMono
        !hamiltonian elemets for Jdown; JdTB stands for Jdown two-body matrix operator
        
        real(kind=4), allocatable, dimension(:,:) :: VHamiltonian !Valence Hamiltonian matrix
    !!!valence energy variables end!!!
    
    real(kind=4), allocatable, dimension(:) :: TEnergy !total energy
    real, allocatable, dimension(:,:) :: TState !total state
    real(kind=4), allocatable, dimension(:,:) :: THamiltonian !whole set of hamiltonian in this model

    !!!M1 transition part begin!!!
        real :: g_proton, g_neutron, g_rotor !g-factors
        real :: g_factor !g_factor=g_proton-g_rotor for proton; g_factor=g_neutron-g_rotor for neutron
        real :: BMRes !final result of B(MI)

        double precision :: MIniState, MFinState !initial I and final I in magnetic transition

        real, allocatable, dimension(:) :: AJplus, AJz, AJminus
        !J matrix of J_{\mu} for part A in M1 transition, \mu=-1, 0, 1

        real, allocatable, dimension(:) :: BJplus, BJz, BJminus
        !J matrix of J_{\mu} for part B in M1 transition, \mu=-1, 0, 1

        real, allocatable, dimension(:) :: CJplus, CJz, CJminus
        !J matrix of J_{\mu} for part C in M1 transition, \mu=-1, 0, 1

        real, allocatable, dimension(:) :: DJplus, DJz, DJminus
        !J matrix of J_{\mu} for part D in M1 transition, \mu=-1, 0, 1
        
        real(kind=4), allocatable, dimension(:,:) :: MIniWF, MFinWF 
        !storing initial and final wave function for transition in magnetic transition
        real, allocatable, dimension(:,:) :: MIniNSCVector, MFinNSCVector, MIniISCVector, MFinISCVector
        !storing inital and final strong coupling basis indexes for K and \phi in magnetic transition
        real, allocatable, dimension(:,:) :: BMResArray 
        !calculated result of B(MI) for different bra&ket in magnetic transition

        double precision :: CGAplus, CGAz, CGAminus
        double precision :: CGBplus, CGBz, CGBminus
        double precision :: CGCplus, CGCz, CGCminus
        double precision :: CGDplus, CGDz, CGDminus
        double precision, dimension(3000,3000,3) :: CGAArray, CGBArray, CGCArray, CGDArray
        !storing CG Factors only for debugging 
    !!!M1 transition part end!!!

    !!!E2 transition part begin!!!
        real :: Q_intr !intrinsic electric quadrupole moment
        real :: BERes

        double precision :: EIniState, EFinState !initial I and final I in electric transition
        double precision :: CG20, CG2P2, CG2N2
        !CG Coefficient result for (IK20|I'K'),(IK22|I'K'),(IK2-2|I'K'), respectively
        double precision :: ICG20, ICG2P2, ICG2N2
        !CG Coefficient result for (IK20|I'-K'),(IK22|I'-K'),(IK2-2|I'-K'), respectively

        real(kind=4), allocatable, dimension(:,:) :: EIniWF, EFinWF
        !storing initial and final wave function for transition in electric transition
        real, allocatable, dimension(:,:) :: EIniNSCVector, EFinNSCVector, EIniISCVector, EFinISCVector
        !storing inital and final strong coupling basis indexes for K and \phi in electric transition
        real, allocatable, dimension(:,:) :: BEResArray 
        !calculated result of B(EI) for different bfa&ket
    !!!E2 transition part end!!!

    !!!angular momentum part begin!!!
        real :: AMState
        !value of I to be calculated for <R^2>
        real :: R1Res, R2Res, R3Res
        !final result of <R_1^2>, <R_1^2>, <R_3^2>
        real :: PJ1Res, PJ2Res, PJ3Res, NJ1Res, NJ2Res, NJ3Res
        !final result of <J_{pk}^2> and <J_{nk}^2>

        real :: RNJUDSMono, RNJUDSPMono, RNJUDSTB
        !normal J matrix for J_+^2+J_-^2
        real :: RNJxySMono, RNJxySPMono, RNJxySTB
        !normal J matirx for J^2-J_3^2
        real :: RNJUp, RNJDown
        !normal J matrix for J_+ and J_-
        real :: RNJ3, RNJ3SMono, RNJ3SPMono, RNJ3STB
        !normal J matrix for J_z and J_z^2

        real :: RIJUDSMono, RIJUDSPMono, RIJUDSTB
        !inverse J matrix for J_+^2+J_-^2
        real :: RIJxySMono, RIJxySPMono, RIJxySTB
        !inverse J matirx for J^2-J_3^2
        real :: RIJUp, RIJDown
        !inverse J matrix for J_+ and J_-
        real :: RIJ3, RIJ3SMono, RIJ3SPMono, RIJ3STB
        !inverse J matrix for J_z and J_z^2

        real :: PNJUDSMono, PNJUDSPMono, PNJUDSTB
        !proton normal J matrix for J_+^2+J_-^2
        real :: PNJxySMono, PNJxySPMono, PNJxySTB
        !proton normal J matrix for J+J- + J-J+
        real :: PNJ3SMono, PNJ3SPMono, PNJ3STB
        !proton normal J matrix for J_3^2

        real :: PIJUDSMono, PIJUDSPMono, PIJUDSTB
        !proton inverse J matrix for J_+^2+J_-^2
        real :: PIJxySMono, PIJxySPMono, PIJxySTB
        !proton inverse J matrix for J+J- + J-J+
        real :: PIJ3SMono, PIJ3SPMono, PIJ3STB
        !proton inverse J matrix for J_3^2

        real :: NNJUDSMono, NNJUDSPMono, NNJUDSTB
        !neutron normal J matrix for J_+^2+J_-^2
        real :: NNJxySMono, NNJxySPMono, NNJxySTB
        !neutron normal J matrix for J+J- + J-J+
        real :: NNJ3SMono, NNJ3SPMono, NNJ3STB
        !neutron normal J matrix for J_3^2

        real :: NIJUDSMono, NIJUDSPMono, NIJUDSTB
        !neutron inverse J matrix for J_+^2+J_-^2
        real :: NIJxySMono, NIJxySPMono, NIJxySTB
        !neutron inverse J matrix for J+J- + J-J+
        real :: NIJ3SMono, NIJ3SPMono, NIJ3STB
        !neutron inverse J matrix for J_3^2

        real, allocatable, dimension(:,:) :: R1Array, R2Array, R3Array
        !array value of <R_1^2>, <R_2^2>, <R_3^2> for each (bra,ket)
        real, allocatable, dimension(:,:) :: PJ1Array, PJ2Array, PJ3Array
        real, allocatable, dimension(:,:) :: NJ1Array, NJ2Array, NJ3Array
        !array value of <J_1^2>, <J_2^2>, <J_3^2> for each (bra,ket) for proton and neutron

        real, allocatable, dimension(:,:) :: AMWF
        !wave function for mean value
        real, allocatable, dimension(:,:) :: AMNSCVector, AMISCVector
    !!!angular momentum part end!!!

    !!!quasi-particle transformation begin!!!
        real, allocatable, dimension(:) :: lambda, Delta
        !Fermi energy and pairing energy gap for each single j shell
        real, allocatable, dimension(:,:) :: QuasiEnergy
        !quasi-particle energy for each \Omega in each single j shell 
        real, allocatable, dimension(:,:) :: QuasiU, QuasiV
        !u_{\mu}, u_{\nu}, v_{\mu} and v_{\nu} from BCS transformation
        !for each \Omega in each single j shell 
    !!!quasi-particle transformation end!!!
   
end module PhysicalQuantities

module DiagonalVariables

    use Constants

    implicit none

    !program process quantities
    integer :: IROT
    real, dimension(10000) :: B, Z
    real, dimension(10000, 10000) :: C !2-dimension matrix C is actually for the EnergyAlignment subroutine

    !judgement values
    integer :: EV = 1
    integer :: NEV = 0

    !physics quantities
    
end module DiagonalVariables

module FunctionParameters

    use Constants

    integer :: PSNumber !particle state serial number
    integer :: FactorialRes !result of a factorial(subroutine)
    integer :: CombinationRes !result of a combination(subroutine)

    double precision :: IniK, FinK

end module FunctionParameters

module CGParameters

    use Constants

    integer :: CombinationRes
    integer :: Factorial123, Factorial11, Factorial22, Factorial321, Factorial312
    !factorial results in CG coefficients
    real :: SumMuRes !sum of every allowable cycle indicator into (-1) formula
    double precision ::  CGResult
    
end module CGParameters


program PRM

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use CGParameters

    implicit none

    !!!for sSYEV begin!!!
        integer(kind=8), parameter :: LWMAX=100000000
        integer(kind=8) :: INFO, LWORK
        real(kind=4), dimension(LWMAX) :: work
    !!!for sSYEV end!!! 

    ! write(*,*) 'total spin'
    ! read(*,*) IStart, IEnd

    ! write(*,*) 'number of single j shell'
    ! read(*,*) tau3_max
    
    ! I=10
    tau3_max=2

    if (.not. allocated(N)) allocate(N(tau3_max))
    if (.not. allocated(MOI)) allocate(MOI(3))
    if (.not. allocated(A)) allocate(A(3))
    if (.not. allocated(Tau3PNum)) allocate(Tau3PNum(10))
    if (.not. allocated(PorN)) allocate(PorN(50))
    if (.not. allocated(ParOrHole)) allocate(ParOrHole(50))
    if (.not. allocated(jsp)) allocate(jsp(tau3_max))
    if (.not. allocated(DeformC)) allocate(DeformC(10))
    if (.not. allocated(countsp_max)) allocate(countsp_max(tau3_max))
    if (.not. allocated(StateSequence)) allocate(StateSequence(tau3_max))
    if (.not. allocated(lambda)) allocate(lambda(tau3_max))
    if (.not. allocated(Delta)) allocate(Delta(tau3_max))

    call Input
    write(*,*) 'input succeeded'

    if (.not. allocated(NilVector)) allocate(NilVector(tau3_max,countsp_max(tcount_max),countsp_max(tcount_max)))
    if (.not. allocated(EAlign)) allocate(EAlign(10,10,countsp_max(tcount_max)))
    if (.not. allocated(CAlign)) allocate(CAlign(10,10,countsp_max(tcount_max),countsp_max(tcount_max)))
    if (.not. allocated(C_nu)) allocate(C_nu(10,-20:20,20))

    do tau3=1, tau3_max
        call SingleParticleEnergy
    end do 
    write(*,*) 'single particle succeeded'

    if (.not. allocated(DecVector)) allocate(DecVector(2,-countsp_max(tcount_max):countsp_max(tcount_max)))
    if (.not. allocated(DecoupleABSum)) allocate(DecoupleABSum(tau3_max,20,20))
    if (.not. allocated(DecoupleCSum)) allocate(DecoupleCSum(tau3_max,20,20))
    if (.not. allocated(DecoupleDSum)) allocate(DecoupleDSum(tau3_max,20,20))
    if (.not. allocated(DecoupleESum)) allocate(DecoupleESum(tau3_max,20,20))
    if (.not. allocated(DecoupleFSum)) allocate(DecoupleFSum(tau3_max,20,20))
    if (.not. allocated(DecoupleGSum)) allocate(DecoupleGSum(tau3_max,20,20))

    call DecouplingFactor
    write(*,*) 'decoupling factor succeeded'

    if (.not. allocated(QuasiEnergy)) allocate(QuasiEnergy(tau3_max,countsp_max(tcount_max)))
    if (.not. allocated(QuasiU)) allocate(QuasiU(tau3_max,countsp_max(tcount_max)))
    if (.not. allocated(QuasiV)) allocate(QuasiV(tau3_max,countsp_max(tcount_max)))

    if (IfQuasi==1) then
        call QuasiTrans
    end if
    write(*,*) 'quasi-particle transformation succeeded'

    if (.not. allocated(SPVector)) allocate(SPVector(10,60,DIMA))
    if (.not. allocated(ISPVector)) allocate(ISPVector(10,60,DIMA))
    if (.not. allocated(PNumN)) allocate(PNumN(10,DIM))
    if (.not. allocated(PNumI)) allocate(PNumI(10,DIM))
    if (.not. allocated(NStateCoupled)) allocate(NStateCoupled(10,DIMA))
    if (.not. allocated(IStateCoupled)) allocate(IStateCoupled(10,DIMA))
    if (.not. allocated(StateNormal)) allocate(StateNormal(tau3_max))
    if (.not. allocated(BSPVector)) allocate(BSPVector(10,60))
    if (.not. allocated(KSPVector)) allocate(KSPVector(10,60))
    if (.not. allocated(JShellIndex)) allocate(JShellIndex(tau3_max,DIM))

    call ParticleSort
    write(*,*) 'particle sort succeeded'

    if (.not. allocated(DimensionArray)) allocate(DimensionArray(50))

    ! !!!calculate energy begin!!! 
    !     Iloop: &
    !     do I=9, 21, 1
    !     countr_max = 2*I+1

    !     call CalDimension
    !     write(*,*) "dimension check=", DimensionNumber
            
    !         if (.not. allocated(KNumber)) allocate(KNumber(DimensionNumber))
    !         if (.not. allocated(NSCVector)) allocate(NSCVector(2,DimensionNumber))
    !         if (.not. allocated(ISCVector)) allocate(ISCVector(2,DimensionNumber))      

    !         call NSCBasisSort
            
    !         write(*,*) 'strong coupling basis succeeded'

    !         if (.not. allocated(IHamiltonian)) allocate(IHamiltonian(DimensionNumber,DimensionNumber))

    !         call IntrinsicEnergy 
    !         write(*,*) 'intrinsic succeeded'

    !         if (.not. allocated(VHamiltonian)) allocate(VHamiltonian(DimensionNumber,DimensionNumber))

    !         call ValenceEnergy
    !         write(*,*) 'valence succeeded'

    !         if (.not. allocated(RHamiltonian)) allocate(RHamiltonian(DimensionNumber,DimensionNumber))
    !         if (.not. allocated(REnergy)) allocate(REnergy(DimensionNumber))
    !         if (.not. allocated(RState)) allocate(RState(DimensionNumber,DimensionNumber))

    !         call RotorEnergy
    !         write(*,*) 'rotor succeeded'

    !         if (.not. allocated(THamiltonian)) allocate(THamiltonian(DimensionNumber,DimensionNumber))
    !         if (.not. allocated(TEnergy)) allocate(TEnergy(DimensionNumber))

    !         THamiltonian = 0
    !         do bra=1, DimensionNumber
    !         do ket=bra, DimensionNumber
    !             THamiltonian(bra,ket) = IHamiltonian(bra,ket) & 
    !             + RHamiltonian(bra,ket) &
    !             + VHamiltonian(bra,ket)
    !         end do
    !         end do

    !         open(unit=Toutput, file='THamiltonian.dat')
    !             ! do bra=1, DimensionNumber
    !             !     write(Toutput,'(100000f12.6)') THamiltonian(bra,:)
    !             ! end do
    !             do bra=1, DimensionNumber
    !             do ket=1, DimensionNumber
    !                 if (THamiltonian(bra,ket)/=0) write(Toutput,*) THamiltonian(bra,ket)
    !             end do
    !             end do
    !         close(unit=Toutput)

    !         ! call Diagonalization(DimensionNumber,DimensionNumber,EV,THamiltonian,TEnergy,TState,IROT,B,Z)

    !         ! call EnergyAlignment(DimensionNumber,DimensionNumber,NEV,TEnergy,C)

    !         INFO=-1
    !         LWORK=-1
    !         call sSYEV('Vectors', 'Upper', DimensionNumber, THamiltonian, DimensionNumber, TEnergy, WORK, LWORK, INFO)
    !         LWORK = min(LWMAX, int(work(1)))
    !         !solve eigen problem!
    !         call sSYEV('Vectors', 'Upper', DimensionNumber, THamiltonian, DimensionNumber, TEnergy, WORK, LWORK, INFO)
    !         !check for convergence!
    !         if (INFO.GT.0) then
    !         write(*,*) 'the alogrithm failed to compute eigenvalues.'
    !         stop
    !         end if

    !         write(*,*) "I=", I
    !         write(*,*) "Total Energy"
    !         write(*,*) TEnergy(1:10)

    !         TFileName = 'WaveFunction/TWaveFunction'//achar(int(I+64))//'.dat'
    !         open(unit=Toutput,file=TFileName)
    !             do bra=1, DimensionNumber
    !                 write(Toutput,'(100000f12.8)') THamiltonian(bra,:)
    !             end do
    !         close(unit=Toutput)

    !         write(*,*) 'totally succeeded' 

    !         if (.not. allocated(KNumber)) allocate(KNumber(DimensionNumber))
    !         if (.not. allocated(NSCVector)) allocate(NSCVector(2,DimensionNumber))
    !         if (.not. allocated(ISCVector)) allocate(ISCVector(2,DimensionNumber)) 

    !         deallocate(TEnergy)
    !         deallocate(THamiltonian)
    !         deallocate(RState)
    !         deallocate(REnergy)
    !         deallocate(RHamiltonian)
    !         deallocate(VHamiltonian)
    !         deallocate(IHamiltonian)
    !         deallocate(ISCVector)
    !         deallocate(NSCVector)
    !         deallocate(KNumber)

    !     end do Iloop

    ! !!!calculate energy end!!!

    ! DimensionArray(12) = 2925
    ! DimensionArray(13) = 3150
    ! DimensionArray(14) = 3375
    ! DimensionArray(15) = 3600
    ! DimensionArray(16) = 3825
    ! DimensionArray(17) = 4050
    ! DimensionArray(18) = 4275
    ! DimensionArray(19) = 4500
    ! DimensionArray(20) = 4725
    ! DimensionArray(21) = 4950
    ! DimensionArray(22) = 5175

    ! DimensionArray(13) = 5544
    ! DimensionArray(14) = 5940
    ! DimensionArray(15) = 6336
    ! DimensionArray(16) = 6732
    ! DimensionArray(17) = 7128
    ! DimensionArray(18) = 7524
    ! DimensionArray(19) = 7920
    ! DimensionArray(20) = 8316
    ! DimensionArray(21) = 8712

    ! DimensionArray(12) = 1521
    ! DimensionArray(13) = 1638
    ! DimensionArray(14) = 1755
    ! DimensionArray(15) = 1872
    ! DimensionArray(16) = 1989
    ! DimensionArray(17) = 2106
    ! DimensionArray(18) = 2223
    ! DimensionArray(19) = 2340
    ! DimensionArray(20) = 2457
    ! DimensionArray(21) = 2574
    ! DimensionArray(22) = 2691

    DimensionArray(9) = 684
    DimensionArray(10) = 756
    DimensionArray(11) = 828
    DimensionArray(12) = 900
    DimensionArray(13) = 972
    DimensionArray(14) = 1044
    DimensionArray(15) = 1116
    DimensionArray(16) = 1188
    DimensionArray(17) = 1260
    DimensionArray(18) = 1332
    DimensionArray(19) = 1404
    DimensionArray(20) = 1476
    DimensionArray(21) = 1548

    ! do sequence=1, TotalStates
    !     write(*,'(100I2)') SPVector(1,1:countsp_max(1),sequence),SPVector(2,1:countsp_max(2),sequence)
    ! end do

    ! stop

    !!!calculate magenetic transitions begin!!!

        MTranloop: &
        do MIniState=10, 21, 1
        MFinState = MIniState-1
        
        allocate(MIniWF(DimensionArray(int(MIniState)),DimensionArray(int(MIniState))))
        allocate(MFinWF(DimensionArray(int(MFinState)),DimensionArray(int(MFinState))))
        allocate(MIniNSCVector(2,DimensionArray(int(MIniState))))
        allocate(MFinNSCVector(2,DimensionArray(int(MFinState))))

        allocate(AJplus(10))
        allocate(AJz(10))
        allocate(AJminus(10))
        allocate(BJplus(10))
        allocate(BJz(10))
        allocate(BJminus(10))
        allocate(CJplus(10))
        allocate(CJz(10))
        allocate(CJminus(10))
        allocate(DJplus(10))
        allocate(DJz(10))
        allocate(DJminus(10))

        allocate(BMResArray(DimensionArray(int(MFinState)),DimensionArray(int(MIniState))))

        call MTransition

        deallocate(BMResArray)

        deallocate(DJminus)
        deallocate(DJz)
        deallocate(DJplus)
        deallocate(CJminus)
        deallocate(CJz)
        deallocate(CJplus)
        deallocate(BJminus)
        deallocate(BJz)
        deallocate(BJplus)
        deallocate(AJminus)
        deallocate(AJz)
        deallocate(AJplus)

        deallocate(MFinNSCVector)
        deallocate(MIniNSCVector)
        deallocate(MFinWF)
        deallocate(MIniWF)

        end do MTranloop

    !!!calculate magnetic transitions end!!!

    ! !!!calculate electric transitions begin!!!

    !     ETranloop: &
    !     do EIniState=15.5, 21.5, 1
    !     EFinState = EIniState-2

    !     ! EIniState = 16.5
    !     ! EFinState = 14.5

    !     allocate(EIniWF(DimensionArray(int(EIniState)),DimensionArray(int(EIniState))))
    !     allocate(EFinWF(DimensionArray(int(EFinState)),DimensionArray(int(EFinState))))
    !     allocate(EIniNSCVector(2,DimensionArray(int(EIniState))))
    !     allocate(EFinNSCVector(2,DimensionArray(int(EFinState))))
    !     allocate(EIniISCVector(2,DimensionArray(int(EIniState))))
    !     allocate(EFinISCVector(2,DimensionArray(int(EFinState))))

    !     allocate(BEResArray(DimensionArray(int(EFinState)),DimensionArray(int(EIniState))))

    !     call ETransition

    !     deallocate(BEResArray)

    !     deallocate(EFinISCVector)
    !     deallocate(EIniISCVector)
    !     deallocate(EFinNSCVector)
    !     deallocate(EIniNSCVector)
    !     deallocate(EFinWF)
    !     deallocate(EIniWF)

    !     end do ETranloop

    ! !!!calculate electric transitions end!!!

    ! !!!calculate angular momentum begin!!!

    !     AngMomloop: &
    !     do AMState=13.5, 22.5, 1

    !     allocate(AMWF(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(AMNSCVector(2,DimensionArray(int(AMState))))
    !     allocate(AMISCVector(2,DimensionArray(int(AMState))))

    !     allocate(R1Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(R2Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(R3Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))

    !     allocate(PJ1Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(PJ2Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(PJ3Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))

    !     allocate(NJ1Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(NJ2Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))
    !     allocate(NJ3Array(DimensionArray(int(AMState)),DimensionArray(int(AMState))))

    !     TFileName='WaveFunction/TWaveFunction'//achar(int(AMState+64))//'.dat'
    !     open(unit=Toutput,file=TFileName)
    !         do bra=1, DimensionArray(int(AMState))
    !             read(Toutput,'(100000f12.8)') AMWF(bra,:)
    !         end do
    !     close(unit=Toutput)

    !     SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(AMState+64))//'.dat'
    !     open(unit=Toutput,file=SCFileName)
    !         do bra=1, DimensionArray(int(AMState))
    !             read(Toutput,*) AMNSCVector(:,bra)
    !         end do
    !     close(unit=Toutput)

    !     SCFileName = 'ISCSerialNum/ISCSerialNum'//achar(int(AMState+64))//'.dat'
    !     open(unit=Toutput,file=SCFileName)
    !         do bra=1, DimensionArray(int(AMState))
    !             read(Toutput,*) AMISCVector(:,bra)
    !         end do
    !     close(unit=Toutput)

    !         ! call AMomentumR

    !         call AMomentumJ

    !     deallocate(NJ3Array)
    !     deallocate(NJ2Array)
    !     deallocate(NJ1Array)

    !     deallocate(PJ3Array)
    !     deallocate(PJ2Array)
    !     deallocate(PJ1Array)

    !     deallocate(R3Array)
    !     deallocate(R2Array)
    !     deallocate(R1Array)

    !     deallocate(AMISCVector)
    !     deallocate(AMNSCVector)
    !     deallocate(AMWF)

    !     end do AngMomloop

    ! !!!calculate angular momentum end!!!

    deallocate(DimensionArray)

    deallocate(JShellIndex)
    deallocate(BSPVector)
    deallocate(KSPVector)
    deallocate(StateNormal)
    deallocate(IStateCoupled)
    deallocate(NStateCoupled)
    deallocate(PNumI)
    deallocate(PNumN)
    deallocate(ISPVector)
    deallocate(SPVector)

    deallocate(QuasiV)
    deallocate(QuasiU)
    deallocate(QuasiEnergy)

    deallocate(DecoupleGSum)
    deallocate(DecoupleFSum)
    deallocate(DecoupleESum)
    deallocate(DecoupleDSum)
    deallocate(DecoupleCSum)
    deallocate(DecoupleABSum)

    deallocate(C_nu)
    deallocate(CAlign)
    deallocate(EAlign)
    deallocate(NilVector)

    deallocate(Delta)
    deallocate(lambda)
    deallocate(StateSequence)
    deallocate(countsp_max)
    deallocate(DeformC)
    deallocate(jsp)
    deallocate(ParOrHole)
    deallocate(PorN)
    deallocate(Tau3PNum)
    deallocate(A)
    deallocate(MOI)
    deallocate(N)

end program PRM