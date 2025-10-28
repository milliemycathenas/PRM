subroutine SingleParticleEnergy

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    !!!for sSYEV begin!!!
        integer(kind=8), parameter :: LWMAX=10000
        integer(kind=8) :: INFO, LWORK
        real(kind=4), dimension(LWMAX) :: work
    !!!for sSYEV end!!!

    SPDimNum = 0.5*countsp_max(tau3)
 
    allocate(SPHamiltonian(SPDimNum,SPDimNum))
    allocate(SPEnergy(SPDimNum))

    NilVectorloop: &
    do count=1, jsp(tau3)+0.5, 1
        NilVector(tau3,1,count)=jsp(tau3)
        if ( mod(jsp(tau3)+0.5, 2.)==0 ) then
            NilVector(tau3,2,count)=-jsp(tau3)+2*(count-1)
        else
            NilVector(tau3,2,count)=-jsp(tau3)+1+2*(count-1)
        end if
    end do NilVectorloop

    SPHamiltonian = 0
    SPEnergy = 0
     
    HspColumnloop:&
    do bra = 1, SPDimNum
        HspRowloop: &
        do ket = bra, SPDimNum

            jsp3_squared = NilVector(tau3,2,ket) * NilVector(tau3,2,ket)
            jsp_squared = NilVector(tau3,1,ket) * (NilVector(tau3,1,ket)+1)
            jspup_squared = sqrt( (NilVector(tau3,1,ket)-NilVector(tau3,2,ket)) * (NilVector(tau3,1,ket)+NilVector(tau3,2,ket)+1) &
                            & * (NilVector(tau3,1,ket)-NilVector(tau3,2,ket)-1) * (NilVector(tau3,1,ket)+NilVector(tau3,2,ket)+2) )
            jspdown_squared = sqrt( (NilVector(tau3,1,ket)+NilVector(tau3,2,ket)) * (NilVector(tau3,1,ket)-NilVector(tau3,2,ket)+1) &
                            & * (NilVector(tau3,1,ket)+NilVector(tau3,2,ket)-1) * (NilVector(tau3,1,ket)-NilVector(tau3,2,ket)+2) )
                                  
            if (NilVector(tau3,2,ket) == NilVector(tau3,2,bra)) then
                SPHamiltonian(bra,ket) = 0.5 * DeformC(tau3) * (cos(gamma_rad) * (jsp3_squared - (jsp_squared/3.)))
            else if ( NilVector(tau3,2,ket) == NilVector(tau3,2,bra) +2 ) then
                SPHamiltonian(bra,ket) = 0.5 * DeformC(tau3) * (sin(gamma_rad) * jspdown_squared)/(2*sqrt(3.))
            else 
                SPHamiltonian(bra,ket) = 0
            end if 
                
            !SPHamiltonian(ket,bra) = SPHamiltonian(bra,ket)

        end do HspRowloop
    end do HspColumnloop

    do ket=2, SPDimNum
    do bra=ket+1, SPDimNum
        SPHamiltonian(bra,ket)=0
    end do
    end do

    !call Diagonalization(DIM, countsp_max(tau3), EV, SPHamiltonian, SPEnergy, SPState, IROT, B, Z)

    INFO=-1
    LWORK=-1
    call sSYEV('Vectors', 'Upper', SPDimNum, SPHamiltonian, SPDimNum, SPEnergy, WORK, LWORK, INFO)
    LWORK = min(LWMAX, int(work(1)))
    !solve eigen problem!
    call sSYEV('Vectors', 'Upper', SPDimNum, SPHamiltonian, SPDimNum, SPEnergy, WORK, LWORK, INFO)
    !check for convergence!
    if (INFO.GT.0) then
    write(*,*) 'the alogrithm failed to compute eigenvalues.'
    stop
    end if

    do count=1, countsp_max(tau3)/2
        EAlign(tau3,1,count)=SPEnergy(count)
    end do

    do row=1, countsp_max(tau3)/2
    do column=1, countsp_max(tau3)/2
        CAlign(tau3,1,row,column)=SPHamiltonian(row,column)
    end do
    end do

    SPEFileName = 'SPEnergy'//achar(tau3+64)//'.dat'
    ! SPSFileName = 'SPState'//achar(tau3+64)//'.dat'

    open(unit=SPoutput, file=SPEFileName)
        SPWriteInloopE: &
        do count=1, SPDimNum
            write(SPoutput, *) SPEnergy(count)
        end do SPWriteInloopE
    close(unit=SPoutput)

    deallocate(SPEnergy)
    deallocate(SPHamiltonian)

    return

end subroutine SingleParticleEnergy



    
