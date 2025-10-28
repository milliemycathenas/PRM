subroutine IntrinsicEnergy
    
    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters
    
    implicit none
 
    IHamiltonian=0
    Hiloop: &
    do bra=1, DimensionNumber
    do ket=1, DimensionNumber

        if (bra==ket) then
    
        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)

                if (IfQuasi==0) then
                    IHamiltonian(bra,ket) = IHamiltonian(bra,ket) &
                    + EAlign(tau3,2,count)*SPVector(tau3,count,NSCVector(2,ket))
                end if

                if (IfQuasi==1) then
                    IHamiltonian(bra,ket) = IHamiltonian(bra,ket) &
                    + QuasiEnergy(tau3,count)*SPVector(tau3,count,NSCVector(2,ket))
                end if
            end do
        end do

        end if

    end do
    end do Hiloop

    open(unit=Ioutput, file='IHamiltonian.dat')
        ! do bra=1, DimensionNumber
        !     write(Ioutput,'(100000f12.6)') IHamiltonian(bra,:)
        ! end do
        do bra=1, DimensionNumber
        do ket=1, DimensionNumber
            if (IHamiltonian(bra,ket)/=0) write(Ioutput,*) IHamiltonian(bra,ket)
        end do
        end do
    close(unit=Ioutput)

end subroutine IntrinsicEnergy
