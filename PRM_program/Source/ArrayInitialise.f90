subroutine ArrayInitialise

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    !!!input parameters begin!!!
        DeformC = 0
        jsp = 0
        countsp_max = 0
        MOI = 0
        A = 0
    !!!input parameters end!!!

    !!!basis vector arrangement begin!!!
        Combination = 0
        
        Tau3PNum = 0
        KNumber = 0

        !ParticleVector = 0
        SPVector = 0
        ISPVector = 0

        NStateCoupled = 0
        IStateCoupled = 0

        PNumN = 0
        PNumI = 0

        NSCVector = 0
        ISCVector = 0
    !!!basis vector arrangement end!!!

    !!!single particle hamiltonian begin!!!
        ParOrHole = 0
        NilVector = 0
        EAlign = 0
        CAlign = 0
        C_nu = 0
        !SPHamiltonian = 0
        !SPEnergy = 0
        QuasiEnergy = 0
        SPState = 0
    !!!single particle hamiltonian end!!!

    
    !!!intrinsic hamiltonian begin!!!
        IHamiltonian = 0
        IEnergy = 0
        IState = 0
    !!!intrinsic hamiltonian end!!!
    
    !!!rotor hamiltonian begin!!!
        RHamiltonian = 0
    !!!rotor hamiltonian end!!!

    !!!decoupling factor begin!!!
        DecVector = 0

        DecoupleAB = 0
        DecoupleC = 0
        DecoupleD = 0
        DecoupleE = 0
        DecoupleF = 0
        DecoupleG = 0

        DecoupleABSum = 0
        DecoupleCSum = 0
        DecoupleDSum = 0
        DecoupleESum = 0
        DecoupleFSum = 0
        DecoupleGSum = 0
    !!!decoupling factor end!!!  

    !!!valence hamiltonian begin!!!
        BSPVector = 0
        KSPVector = 0

        VHamiltonian = 0
    !!!valence hamiltonian end!!!

    !!!total hamiltonian begin!!!
        THamiltonian = 0
        TEnergy = 0
        TState = 0
    !!!total hamiltonian end!!!

    !!!diagonalization parameter begin!!!
        B = 0
        Z = 0
        C = 0    
    !!!diagonalization parameter end!!!

end subroutine ArrayInitialise