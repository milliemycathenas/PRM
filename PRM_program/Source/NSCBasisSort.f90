subroutine NSCBasisSort
!occupied particle number

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables

    implicit none 

    if (IfTrunc==1) then
        E_min = 0
        do tau3=1, tau3_max
            !normal state
            do count=1, int(Tau3PNum(tau3)/2)
                E_min = E_min + EAlign(tau3,2,count)
            end do
            !inverse state
            do count=1, Tau3PNum(tau3)-int(Tau3PNum(tau3)/2)
                E_min = E_min + EAlign(tau3,2,count)
            end do
        end do
    end if
    
    bra = 0
    KNumber = 0
    do sequence=1, TotalStates
        KNumber(sequence) = 0
        do count=1, countr_max
            SignPass = -I + (count-1)

            PassMark=0
            do tau3=1, tau3_max
                PassMark = PassMark + (PNumN(tau3,sequence) - PNumI(tau3,sequence))
            end do

            if (IfTrunc==1) then
            
                E_difference = 0
                E_MPC = 0
                do tau3=1, tau3_max
                    do serialnum=1, countsp_max(tau3)
                        E_MPC = E_MPC + EAlign(tau3,2,serialnum) &
                        *SPVector(tau3,serialnum,sequence)
                    end do
                end do
                E_difference = E_MPC - E_min

                D2Symmetry = SignPass - 0.5*PassMark
                if ( mod(D2Symmetry, 2)==0 .and. E_difference<=E_cut) then
                    bra = bra + 1  
                    KNumber(sequence) = KNumber(sequence) + 1 

                    NSCVector(1,bra) = SignPass
                    NSCVector(2,bra) = sequence
                    ISCVector(1,bra) = -SignPass
                    ISCVector(2,bra) = sequence
                end if
            end if

            if (IfTrunc==0) then
                D2Symmetry = SignPass - 0.5*PassMark
                if ( mod(D2Symmetry, 2)==0 ) then
                    bra = bra + 1  
                    KNumber(sequence) = KNumber(sequence) + 1 

                    NSCVector(1,bra) = SignPass
                    NSCVector(2,bra) = sequence
                    ISCVector(1,bra) = -SignPass
                    ISCVector(2,bra) = sequence
                end if
            end if
        end do
    end do 

    DimensionNumber=sum(KNumber)
    write(*,*) 'DIM_trunc=', DimensionNumber

    DimensionArray(int(I)) = DimensionNumber

    ! do bra=1, DimensionNumber
    !     write(*,*) NSCVector(1,bra)
    !     write(*,'(100I2)') SPVector(1,1:countsp_max(1),NSCVector(2,bra)),SPVector(2,1:countsp_max(2),NSCVector(2,bra))
    !     write(*,'(100I2)') SPVector(2,1:countsp_max(2),NSCVector(2,bra))
    ! end do

    ! stop

    SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(I+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionNumber
            write(Toutput,'(2f12.2)') NSCVector(:,bra)
        end do
    close(unit=Toutput)

    SCFileName = 'ISCSerialNum/ISCSerialNum'//achar(int(I+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionNumber
            write(Toutput,'(2f12.2)') ISCVector(:,bra)
        end do
    close(unit=Toutput)
    
end subroutine NSCBasisSort