subroutine RotorEnergy

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    ! !!!for sSYEV begin!!!
    !     integer(kind=8), parameter :: LWMAX=100000000
    !     integer(kind=8) :: INFO, LWORK
    !     real(kind=4), dimension(LWMAX) :: work
    ! !!!for sSYEV end!!! 

    RHamiltonian = 0

    HrColumnloop: &
    do bra = 1, DimensionNumber
        HrRowloop: &
        do ket = bra, DimensionNumber

            ISCVector(1,ket) = -NSCVector(1,ket)

            if( NSCVector(2,ket) == NSCVector(2,bra) ) then

                Irotor3 = NSCVector(1,ket)
                Irotor_squared = I*(I+1)

                Irotorup_squared = sqrt( ( ( I*(I+1) ) - ( NSCVector(1,ket)*(NSCVector(1,ket)+1) ) ) &
                * ( ( I*(I+1) ) - ( (NSCVector(1,ket)+1)*(NSCVector(1,ket)+2) ) ) )

                Irotordown_squared = sqrt( ( ( I*(I+1) ) - ( NSCVector(1,ket)*(NSCVector(1,ket)-1) ) ) &
                * ( ( I*(I+1) ) - ( (NSCVector(1,ket)-1) * (NSCVector(1,ket)-2) ) ) )

                if ( NSCVector(1,ket) == NSCVector(1,bra)+2 ) then
                    RHamiltonian(bra,ket) = 0.25 * A(1) * Irotordown_squared 
                end if
                if ( NSCVector(1,ket) == NSCVector(1,bra)-2 ) then
                    RHamiltonian(bra,ket) = 0.25 * A(1) * Irotorup_squared
                end if
                if ( NSCVector(1,ket) == NSCVector(1,bra) ) then
                    RHamiltonian(bra,ket) = 0.5 * A(2) * (Irotor_squared - Irotor3*Irotor3) + A(3) * (Irotor3*Irotor3)
                end if

            end if

        end do HrRowloop
    end do HrColumnloop

    ! call Diagonalization(DimensionNumber, DimensionNumber, EV, RHamiltonian, REnergy, RState, IROT, B, Z)

    ! INFO=-1
    ! LWORK=-1
    ! call sSYEV('Vectors', 'Upper', DimensionNumber, RHamiltonian, DimensionNumber, REnergy, WORK, LWORK, INFO)
    ! LWORK = min(LWMAX, int(work(1)))
    ! !solve eigen problem!
    ! call sSYEV('Vectors', 'Upper', DimensionNumber, RHamiltonian, DimensionNumber, REnergy, WORK, LWORK, INFO)
    ! !check for convergence!
    ! if (INFO.GT.0) then
    ! write(*,*) 'the alogrithm failed to compute eigenvalues.'
    ! stop
    ! end if

    ! RFileName='REnergy'//achar(int(I)+64)//'.dat'

    ! open(unit=Routput, file=RFileName)
    !     do bra=1, DimensionNumber
    !         write(Routput,*) REnergy(bra) 
    !     end do
    ! close(unit=Routput)

    ! SignPass=0
    ! do ket=1, DimensionNumber-1
    ! do bra=ket+1, DimensionNumber
    !     RHamiltonian(bra,ket)=0
    !     SignPass=SignPass+1
    ! end do
    ! end do

    ! write(*,*) 'sign', SignPass

    open(unit=Routput, file='RHamiltonian.dat')
        ! do bra=1, DimensionNumber
        !     write(Routput,'(100000f12.6)') RHamiltonian(bra,:)
        ! end do
        do bra=1, DimensionNumber
        do ket=1, DimensionNumber
            if (RHamiltonian(bra,ket)/=0) write(Routput,*) RHamiltonian(bra,ket)
        end do
        end do
    close(unit=Routput)

end subroutine RotorEnergy