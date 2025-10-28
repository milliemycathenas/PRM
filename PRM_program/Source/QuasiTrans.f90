subroutine QuasiTrans

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters

    implicit none

    do tau3=1, tau3_max
        do count=1, countsp_max(tau3)
            QuasiEnergy(tau3,count)=sqrt((EAlign(tau3,2,count)-lambda(tau3))**2+Delta(tau3)**2)
        end do
    end do

    do tau3=1, tau3_max
        do count=1, countsp_max(tau3)
            QuasiV(tau3,count)=sqrt(0.5*(1-((EAlign(tau3,2,count)-lambda(tau3))/QuasiEnergy(tau3,count))))
        end do
    end do

    do tau3=1, tau3_max
        do count=1, countsp_max(tau3)
            QuasiU(tau3,count)=sqrt(1-QuasiV(tau3,count)**2)
        end do
    end do
    
    do tau3=1, tau3_max
        do CountMu=1, countsp_max(tau3)
        do CountNu=1, countsp_max(tau3)
            DecoupleABSum(tau3,CountMu,CountNu) = DecoupleABSum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))

            DecoupleCSum(tau3,CountMu,CountNu) = DecoupleCSum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))

            DecoupleDSum(tau3,CountMu,CountNu) = DecoupleDSum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))

            DecoupleESum(tau3,CountMu,CountNu) = DecoupleESum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))

            DecoupleFSum(tau3,CountMu,CountNu) = DecoupleFSum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))

            DecoupleGSum(tau3,CountMu,CountNu) = DecoupleGSum(tau3,CountMu,CountNu) &
            *(QuasiU(tau3,CountMu)*QuasiU(tau3,CountNu)+QuasiV(tau3,CountMu)*QuasiV(tau3,CountNu))
        end do
        end do
    end do

    ! write(*,*) 'after BCS'

    ! do tau3=1, tau3_max
    ! write(*,*) 'tau3', tau3, 'F'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleFSum(tau3,column,1:countsp_max(tau3))
    ! end do
    ! end do


end subroutine QuasiTrans