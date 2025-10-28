subroutine Input

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    ! !!!input parameters begin!!!
    !     write(*,*) 'nuclide mass:'
    !     read(*,*) Mass

    !     hbaromega_0 = 41.*(Mass**(-1/3.))

    !     write(*,*) 'proton number:'
    !     read(*,*) NucProtonNum

    !     write(*,*) 'beta:'
    !     read(*,*) beta

    !     write(*,*) 'gamma:'
    !     read(*,*) gamma
    !     gamma_rad = gamma*pi/180.

    !     do tau3=1, tau3_max
    !         write(*,*) 'Nilsson quantum number N:', tau3
    !         read(*,*) N(tau3)
    !     end do

    !     write(*,*) 'motion of inertia:'
    !     read(*,*) MOI_0

    !     MOIloop: &
    !     do count=1, 3
    !         MOI(count) = MOI_0 * ((sin(gamma_rad - (2*pi*count)/3))**2)
    !     end do MOIloop

    !     A(1) = 1/(2*MOI(1)) - 1/(2*MOI(2))
    !     A(2) = 1/(2*MOI(1)) + 1/(2*MOI(2))
    !     A(3) = 1/(2*MOI(3))

    !     do tau3=1, tau3_max
    !         write(*,*) 'number of particle or hole on one single j shell:', tau3
    !         read(*,*) Tau3PNum(tau3)
    !         write(*,*) 'is this proton or neutron (proton=1, neutron=0):'
    !         read(*,*) PorN(tau3)
    !         write(*,*) 'is this particle or hole (particle=1, hole=-1):'
    !         read(*,*) ParOrHole(tau3)
    !         write(*,*) 'single j shell for this particle:'
    !         read(*,*) jsp(tau3)
    !     end do

    !     do tau3=1, tau3_max
    !         DeformC(tau3) = (sqrt(5./(4.*pi))*1.5*beta)*(1/(jsp(tau3)*(jsp(tau3)+1))) &
    !         *hbaromega_0*(N(tau3)+1.5)
    !         ! DeformC(tau3) = (38.8*(N+1.5)*(Mass**(-1/3.))*beta)/(jsp(tau3)*(jsp(tau3)+1))
    !         if (ParOrHole(tau3)==-1) DeformC(tau3)=-DeformC(tau3)
    !     end do
        
    !     do tau3=1, tau3_max
    !         countsp_max(tau3) = 2*jsp(tau3)+1
    !     end do

    !     tcount_max=1
    !     do tau3=1, tau3_max
    !         if (countsp_max(tau3+1)>countsp_max(tau3)) tcount_max=tau3+1
    !     end do
        
    !     do tau3=1, tau3_max
    !         StateSequence(tau3) = 2**countsp_max(tau3)
    !     end do

    !     write(*,*) 'calculate eletromagnetic transitions? (yes=1, no=0): '
    !     read(*,*) IfEMTrans

    !     if (IfEMTrans==1) then
    !         NucRadius = 1.2*(Mass**(1/3.))
    !         Q_intr = 0.03*(NucRadius**2)*NucProtonNum*beta/sqrt(5*pi)

    !         write(*,*) 'proton g factor:'
    !         read(*,*) g_factor
    !         write(*,*) 'neutron g factor'
    !         read(*,*) g_neutron
    !         write(*,*) 'rotor g factor'
    !         read(*,*) g_rotor
    !     end if
        
    !     write(*,*) 'BCS tranformation? (yes=1, no=0):'
    !     read(*,*) IfQuasi

    !     if (IfQuasi==1) then
    !         do tau3=1, tau3_max
    !             write(*,*) 'Fermi energy for particle on shell', tau3
    !             read(*,*) lambda(tau3)
    !             write(*,*) 'pairing energy gap for particle on shell', tau3
    !             read(*,*) Delta(tau3)
    !         end do
    !     end if

    !     write(*,*) 'if truncate or not? (yes=1, no=0)'
    !     read(*,*) IfTrunc

    !     if (IfTrunc==1) then
    !         write(*,*) 'truncation energy'
    !         read(*,*) E_cut
    !     end if

    ! !!!input parameters end!!!

    ! !!!135Nd begin!!!
    !     Nuclide='Nd'
    !     Mass=135
    !     ! write(NuclideMass,*) Mass
    !     ! write(*,*) NuclideMass

    !     NucProtonNum=60
    !     beta=0.235
    !     gamma=22.4
    !     N(1)=5
    !     N(2)=5
    !     MOI_0=29
    !     Tau3PNum(1)=2
    !     Tau3PNum(2)=1
    !     PorN(1)=1
    !     PorN(2)=0
    !     ParOrHole(1)=1
    !     ParOrHole(2)=-1
    !     jsp(1)=5.5
    !     jsp(2)=5.5

    !     IfTrunc=1
    !     IfEMTrans=1
    !     IfQuasi=0
    !     CoriolisF=1
    !     E_cut=6

    !     hbaromega_0 = 41.*(Mass**(-1/3.))
    !     gamma_rad = gamma*pi/180.

    !     MOIloop: &
    !     do count=1, 3
    !         MOI(count) = MOI_0 * ((sin(gamma_rad - (2*pi*count)/3))**2)
    !     end do MOIloop

    !     A(1) = 1/(2*MOI(1)) - 1/(2*MOI(2))
    !     A(2) = 1/(2*MOI(1)) + 1/(2*MOI(2))
    !     A(3) = 1/(2*MOI(3))

    !     do tau3=1, tau3_max
    !         DeformC(tau3) = (sqrt(5./(4.*pi))*1.5*beta)*(1/(jsp(tau3)*(jsp(tau3)+1))) &
    !         *hbaromega_0*(N(tau3)+1.5)
    !         ! DeformC(tau3) = (38.8*(N+1.5)*(Mass**(-1/3.))*beta)/(jsp(tau3)*(jsp(tau3)+1))
    !         if (ParOrHole(tau3)==-1) DeformC(tau3)=-DeformC(tau3)
    !     end do
        
    !     do tau3=1, tau3_max
    !         countsp_max(tau3) = 2*jsp(tau3)+1
    !     end do

    !     tcount_max=1
    !     do tau3=1, tau3_max
    !         if (countsp_max(tau3+1)>countsp_max(tau3)) tcount_max=tau3+1
    !     end do
        
    !     do tau3=1, tau3_max
    !         StateSequence(tau3) = 2**countsp_max(tau3)
    !     end do

    !     g_rotor=0.44
    !     g_proton=1.11
    !     g_neutron=-0.21
        
    
    !     if (IfEMTrans==1) then
    !         NucRadius = 1.2*(Mass**(1/3.))
    !         Q_intr = 0.03*(NucRadius**2)*NucProtonNum*beta/sqrt(5*pi)
    !     end if
    ! !!!135Nd end!!!

    !!!126Cs begin!!!
        Mass=126
        NucProtonNum=55
        beta=0.26
        gamma=24
        N(1)=5
        N(2)=5
        MOI_0=20
        Tau3PNum(1)=1
        Tau3PNum(2)=1
        PorN(1)=1
        PorN(2)=0
        ParOrHole(1)=1
        ParOrHole(2)=1
        jsp(1)=5.5
        jsp(2)=5.5

        IfTrunc=0
        IfEMTrans=1
        IfQuasi=1
        CoriolisF=1

        lambda(1)=-2.293
        lambda(2)=0.8

        Delta(1)=1
        Delta(2)=1

        hbaromega_0 = 41.*(Mass**(-1/3.))
        gamma_rad = gamma*pi/180.

        MOIloop: &
        do count=1, 3
            MOI(count) = MOI_0 * ((sin(gamma_rad - (2*pi*count)/3))**2)
        end do MOIloop

        A(1) = 1/(2*MOI(1)) - 1/(2*MOI(2))
        A(2) = 1/(2*MOI(1)) + 1/(2*MOI(2))
        A(3) = 1/(2*MOI(3))

        do tau3=1, tau3_max
            DeformC(tau3) = (sqrt(5./(4.*pi))*1.5*beta)*(1/(jsp(tau3)*(jsp(tau3)+1))) &
            *hbaromega_0*(N(tau3)+1.5)
            ! DeformC(tau3) = (38.8*(N+1.5)*(Mass**(-1/3.))*beta)/(jsp(tau3)*(jsp(tau3)+1))
            if (ParOrHole(tau3)==-1) DeformC(tau3)=-DeformC(tau3)
        end do
        
        do tau3=1, tau3_max
            countsp_max(tau3) = 2*jsp(tau3)+1
        end do

        tcount_max=1
        do tau3=1, tau3_max
            if (countsp_max(tau3+1)>countsp_max(tau3)) tcount_max=tau3+1
        end do
        
        do tau3=1, tau3_max
            StateSequence(tau3) = 2**countsp_max(tau3)
        end do

        DeformC(1)=0.3
        DeformC(2)=0.3

        g_rotor=0.44
        g_proton=1.21
        g_neutron=-0.21
        Q_intr=3.5
    !!!126Cs end!!!

    ! !!!124Cs begin!!!
    !     Mass=124
    !     NucProtonNum=55
    !     beta=0.27
    !     gamma=22.1
    !     N(1)=5
    !     N(2)=5
    !     MOI_0=18
    !     Tau3PNum(1)=1
    !     Tau3PNum(2)=1
    !     PorN(1)=1
    !     PorN(2)=0
    !     ParOrHole(1)=1
    !     ParOrHole(2)=1
    !     jsp(1)=5.5
    !     jsp(2)=5.5

    !     IfTrunc=0
    !     IfEMTrans=1
    !     IfQuasi=1

    !     lambda(1)=-2.293
    !     lambda(2)=0.8

    !     Delta(1)=1
    !     Delta(2)=1

    !     hbaromega_0 = 41.*(Mass**(-1/3.))
    !     gamma_rad = gamma*pi/180.

    !     MOIloop: &
    !     do count=1, 3
    !         MOI(count) = MOI_0 * ((sin(gamma_rad - (2*pi*count)/3))**2)
    !     end do MOIloop

    !     A(1) = 1/(2*MOI(1)) - 1/(2*MOI(2))
    !     A(2) = 1/(2*MOI(1)) + 1/(2*MOI(2))
    !     A(3) = 1/(2*MOI(3))

    !     do tau3=1, tau3_max
    !         DeformC(tau3) = (sqrt(5./(4.*pi))*1.5*beta)*(1/(jsp(tau3)*(jsp(tau3)+1))) &
    !         *hbaromega_0*(N(tau3)+1.5)
    !         ! DeformC(tau3) = (38.8*(N+1.5)*(Mass**(-1/3.))*beta)/(jsp(tau3)*(jsp(tau3)+1))
    !         if (ParOrHole(tau3)==-1) DeformC(tau3)=-DeformC(tau3)
    !     end do
        
    !     do tau3=1, tau3_max
    !         countsp_max(tau3) = 2*jsp(tau3)+1
    !     end do

    !     tcount_max=1
    !     do tau3=1, tau3_max
    !         if (countsp_max(tau3+1)>countsp_max(tau3)) tcount_max=tau3+1
    !     end do
        
    !     do tau3=1, tau3_max
    !         StateSequence(tau3) = 2**countsp_max(tau3)
    !     end do

    !     DeformC(1)=0.385
    !     DeformC(2)=0.385

    !     g_rotor=0.44
    !     g_proton=1.21
    !     g_neutron=-0.21
    !     Q_intr=4

    ! !!!124Cs end!!!

    ! !!!105Ag begin!!!
    !     Mass=105
    !     NucProtonNum=47
    !     beta=0.23
    !     gamma=30
    !     N(1)=4
    !     N(2)=5
    !     N(3)=4
    !     MOI_0=32
    !     Tau3PNum(1)=1
    !     Tau3PNum(2)=1
    !     Tau3PNum(3)=1
    !     PorN(1)=1
    !     PorN(2)=0
    !     PorN(3)=0
    !     ParOrHole(1)=-1
    !     ParOrHole(2)=1
    !     ParOrHole(3)=-1
    !     jsp(1)=4.5
    !     jsp(2)=5.5
    !     jsp(3)=2.5

    !     IfTrunc=0
    !     IfEMTrans=0
    !     IfQuasi=0
    !     CoriolisF=0.85
    !     E_cut=6

    !     hbaromega_0 = 41.*(Mass**(-1/3.))
    !     gamma_rad = gamma*pi/180.

    !     MOIloop: &
    !     do count=1, 3
    !         MOI(count) = MOI_0 * ((sin(gamma_rad - (2*pi*count)/3))**2)
    !     end do MOIloop

    !     A(1) = 1/(2*MOI(1)) - 1/(2*MOI(2))
    !     A(2) = 1/(2*MOI(1)) + 1/(2*MOI(2))
    !     A(3) = 1/(2*MOI(3))

    !     do tau3=1, tau3_max
    !         DeformC(tau3) = (sqrt(5./(4.*pi))*1.5*beta)*(1/(jsp(tau3)*(jsp(tau3)+1))) &
    !         *hbaromega_0*(N(tau3)+1.5)
    !         ! DeformC(tau3) = (38.8*(N+1.5)*(Mass**(-1/3.))*beta)/(jsp(tau3)*(jsp(tau3)+1))
    !         if (ParOrHole(tau3)==-1) DeformC(tau3)=-DeformC(tau3)
    !     end do
        
    !     do tau3=1, tau3_max
    !         countsp_max(tau3) = 2*jsp(tau3)+1
    !     end do

    !     tcount_max=1
    !     do tau3=1, tau3_max
    !         if (countsp_max(tau3+1)>countsp_max(tau3)) tcount_max=tau3+1
    !     end do
        
    !     do tau3=1, tau3_max
    !         StateSequence(tau3) = 2**countsp_max(tau3)
    !     end do   
    ! !!!105Ag end!!!

    ! do tau3=1, tau3_max
    !     write(*,*) 'deform C for shell', tau3, ':'
    !     read(*,*) DeformC(tau3)
    ! end do

end subroutine Input
