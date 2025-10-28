subroutine DecouplingFactor

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables

    implicit none

    call ReAlign

    DecoupleABSum=0
    DecoupleCSum=0
    DecoupleDSum=0
    DecoupleESum=0
    DecoupleFSum=0
    DecoupleGSum=0

    if (.not. allocated(DecoupleAB)) allocate(DecoupleAB(tau3_max,-20:20))
    if (.not. allocated(DecoupleC)) allocate(DecoupleC(tau3_max,-20:20))
    if (.not. allocated(DecoupleD)) allocate(DecoupleD(tau3_max,-20:20))
    if (.not. allocated(DecoupleE)) allocate(DecoupleE(tau3_max,-20:20))
    if (.not. allocated(DecoupleF)) allocate(DecoupleF(tau3_max,-20:20))
    if (.not. allocated(DecoupleG)) allocate(DecoupleG(tau3_max,-20:20))

    tau3loop: &
    do tau3=1, tau3_max

    DecVectorloop: &
    do count=1,countsp_max(tau3)
        DecVector(1,count) = jsp(tau3)
        DecVector(2,count) = -jsp(tau3)+(count-1)
    end do DecVectorloop

    DFloop: &
    do CountMu=1, countsp_max(tau3)
        do CountNu=1, countsp_max(tau3)

            DecoupleABloop: &
            do count=1, countsp_max(tau3)
                DecoupleAB(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count-1,CountNu) &
                *sqrt((DecVector(1,count)*(DecVector(1,count)+1)) - (DecVector(2,count)*(DecVector(2,count)-1)))
            end do DecoupleABloop

            DecoupleCloop: &
            do count=1, countsp_max(tau3)
                DecoupleC(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count+1,CountNu) &
                *sqrt((DecVector(1,count)*(DecVector(1,count)+1)) - (DecVector(2,count)*(DecVector(2,count)+1)))
            end do DecoupleCloop

            DecoupleDloop: &
            do count=1, countsp_max(tau3)
                DecoupleD(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count,CountNu)*DecVector(2,count)
            end do DecoupleDloop

            DecoupleEloop: &
            do count=1, countsp_max(tau3)
                DecoupleE(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count,CountNu)*(DecVector(2,count)**2)
            end do DecoupleEloop

            DecoupleFloop: &
            do count=1, countsp_max(tau3)
                DecoupleF(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count,CountNu) &
                *(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)**2))
            end do DecoupleFloop

            DecoupleGloop: &
            do count=1, countsp_max(tau3)
                if (DecVector(1,count)>abs(DecVector(2,count))) then
                    DecoupleG(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count-2,CountNu) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)-2)*(DecVector(2,count)-1)) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)-1)*DecVector(2,count)) &
                    + C_nu(tau3,count,CountMu)*C_nu(tau3,count+2,CountNu) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)+2)*(DecVector(2,count)+1)) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)+1)*DecVector(2,count))
                else if (DecVector(2,count)==DecVector(1,count)) then
                    DecoupleG(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count-2,CountNu) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)-2)*(DecVector(2,count)-1)) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)-1)*DecVector(2,count))
                else if (DecVector(2,count)==-DecVector(1,count)) then
                    DecoupleG(tau3,count) = C_nu(tau3,count,CountMu)*C_nu(tau3,count+2,CountNu) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)+2)*(DecVector(2,count)+1)) &
                    *sqrt(DecVector(1,count)*(DecVector(1,count)+1) - (DecVector(2,count)+1)*DecVector(2,count))
                end if
            end do DecoupleGloop

            DecoupleABSum(tau3,CountMu,CountNu) = sum(DecoupleAB(tau3,1:countsp_max(tau3)))
            DecoupleCSum(tau3,CountMu,CountNu) = sum(DecoupleC(tau3,1:countsp_max(tau3)))
            DecoupleDSum(tau3,CountMu,CountNu) = sum(DecoupleD(tau3,1:countsp_max(tau3)))
            DecoupleESum(tau3,CountMu,CountNu) = sum(DecoupleE(tau3,1:countsp_max(tau3)))           
            DecoupleFSum(tau3,CountMu,CountNu) = sum(DecoupleF(tau3,1:countsp_max(tau3)))
            DecoupleGSum(tau3,CountMu,CountNu) = sum(DecoupleG(tau3,1:countsp_max(tau3)))

        end do
    end do DFloop

    ! write(*,*) 'tau3', tau3, 'AB'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleABSum(tau3,column,1:countsp_max(tau3))
    ! end do

    ! write(*,*) 'tau3', tau3, 'C'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleCSum(tau3,column,1:countsp_max(tau3))
    ! end do

    ! write(*,*) 'tau3', tau3, 'D'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleDSum(tau3,column,1:countsp_max(tau3))
    ! end do

    ! write(*,*) 'tau3', tau3, 'E'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleESum(tau3,column,1:countsp_max(tau3))
    ! end do

    ! write(*,*) 'tau3', tau3, 'F'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleFSum(tau3,column,1:countsp_max(tau3))
    ! end do

    ! write(*,*) 'tau3', tau3, 'G'
    ! do column=1, countsp_max(tau3)
    !     write(*,'(100f12.6)') DecoupleGSum(tau3,column,1:countsp_max(tau3))
    ! end do

    end do tau3loop

    ! stop

    deallocate(DecoupleG)
    deallocate(DecoupleF)
    deallocate(DecoupleE)
    deallocate(DecoupleD)
    deallocate(DecoupleC)
    deallocate(DecoupleAB)

end subroutine DecouplingFactor