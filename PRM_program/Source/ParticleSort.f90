subroutine ParticleSort

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters

    implicit none

    integer :: JShell1, JShell2, JShell3, JShell4, JShell5
    !loop variables for particle organization
    integer, dimension(10,DIM) :: NIDiffer !PNumN(tau3)-PNumI(tau3)

    if (.not. allocated(StateRemained)) allocate(StateRemained(DIM))
    if (.not. allocated(PhiVector)) allocate(PhiVector(100,DIM))
    if (.not. allocated(IPhiVector)) allocate(IPhiVector(100,DIM))
    if (.not. allocated(ParticleVector)) allocate(ParticleVector(DIM))

    do tau3=1, tau3_max
        ParticleVector=0
        StateRemained(tau3)=0

        SortOccupationloop: &
        do sequence=1, StateSequence(tau3)
            serialnum=sequence-1

            BinaryTransformloop: &
            do count=countsp_max(tau3), 1, -1
                ParticleVector(count) = mod(serialnum,2)
                serialnum=serialnum/2
            end do BinaryTransformloop

            if (sum(ParticleVector) == Tau3PNum(tau3)) then
                StateRemained(tau3)=StateRemained(tau3)+1

                do count=1, countsp_max(tau3)
                    SPVector(tau3,count,StateRemained(tau3))=ParticleVector(count)
                end do
            end if
        end do SortOccupationloop
    end do 

    do tau3=1, tau3_max
    do sequence=1, StateRemained(tau3)
        SignPass=countsp_max(tau3)/2
        PNumN(tau3,sequence)=sum(SPVector(tau3,1:SignPass,sequence))
        PNumI(tau3,sequence)=sum(SPVector(tau3,SignPass+1:countsp_max(tau3),sequence))
    end do
    end do

    if (tau3_max==2) then
        PassMark=0
        do column=1, StateRemained(1)
        do row=1, StateRemained(2)
            PassMark=PassMark+1
            JShellIndex(1,PassMark)=column
            JShellIndex(2,PassMark)=row
        end do  
        end do

        TotalStates=0
        do sequence=1, PassMark
            if (PNumN(1,JShellIndex(1,sequence))-PNumI(1,JShellIndex(1,sequence))+ &
            PNumN(2,JShellIndex(2,sequence))-PNumI(2,JShellIndex(2,sequence))>0 &
            .or. (PNumN(1,JShellIndex(1,sequence))-PNumI(1,JShellIndex(1,sequence))+ &
            PNumN(2,JShellIndex(2,sequence))-PNumI(2,JShellIndex(2,sequence))==0  &
            .and. PNumN(1,JShellIndex(1,sequence))-PNumI(1,JShellIndex(1,sequence))>=0)) then
                
                TotalStates=TotalStates+1

                EndPoint=sum(countsp_max(1:1))
                BeginPoint=EndPoint-countsp_max(1)
                do count=1, countsp_max(1)
                    PhiVector(count+BeginPoint,TotalStates)= &
                    SPVector(1,count,JShellIndex(1,sequence))
                end do

                EndPoint=sum(countsp_max(1:2))
                BeginPoint=EndPoint-countsp_max(2)
                do count=1, countsp_max(2)
                    PhiVector(count+BeginPoint,TotalStates)= &
                    SPVector(2,count,JShellIndex(2,sequence))
                end do
            
            end if
        end do
    end if

    if (tau3_max==3) then
        PassMark=0
        do JShell1=1, StateRemained(1)
        do JShell2=1, StateRemained(2)
        do JShell3=1, StateRemained(3)
            PassMark=PassMark+1
            JShellIndex(1,PassMark)=JShell1
            JShellIndex(2,PassMark)=JShell2
            JShellIndex(3,PassMark)=JShell3
        end do
        end do
        end do

        do sequence=1, PassMark
        do tau3=1, tau3_max
            NIDiffer(tau3,sequence)=PNumN(tau3,JShellIndex(tau3,sequence)) &
                                    -PNumI(tau3,JShellIndex(tau3,sequence))
        end do
        end do

        TotalStates=0
        do sequence=1, PassMark
            if (sum(NIDiffer(:,sequence))>0 &
                .or. (sum(NIDiffer(:,sequence))==0 .and. NIDiffer(1,sequence)>0 .and. NIDiffer(2,sequence)<0) &
                .or. (sum(NIDiffer(:,sequence))==0 .and. NIDiffer(1,sequence)>0 .and. NIDiffer(3,sequence)<0) &
                .or. (sum(NIDiffer(:,sequence))==0 .and. NIDiffer(2,sequence)>0 .and. NIDiffer(3,sequence)<0) &
                .or. (sum(NIDiffer(:,sequence))==0 &
                .and. NIDiffer(1,sequence)==0 .and. NIDiffer(2,sequence)==0 .and. NIDiffer(3,sequence)==0)) then

                TotalStates=TotalStates+1

                EndPoint=sum(countsp_max(1:1))
                BeginPoint=EndPoint-countsp_max(1)
                do count=1, countsp_max(1)
                    PhiVector(count+BeginPoint,TotalStates)= &
                    SPVector(1,count,JShellIndex(1,sequence))
                end do

                EndPoint=sum(countsp_max(1:2))
                BeginPoint=EndPoint-countsp_max(2)
                do count=1, countsp_max(2)
                    PhiVector(count+BeginPoint,TotalStates)= &
                    SPVector(2,count,JShellIndex(2,sequence))
                end do

                EndPoint=sum(countsp_max(1:3))
                BeginPoint=EndPoint-countsp_max(3)
                do count=1, countsp_max(3)
                    PhiVector(count+BeginPoint,TotalStates)= &
                    SPVector(3,count,JShellIndex(3,sequence))
                end do
            end if
        end do
    end if

    SPVector=0

    do sequence=1, TotalStates
        do tau3=1, tau3_max
        EndPoint=sum(countsp_max(1:tau3-1))
            do count=1, countsp_max(tau3)
                SPVector(tau3,count,sequence)=PhiVector(count+EndPoint,sequence)
            end do
        end do
    end do

    do sequence=1, TotalStates
    do tau3=1, tau3_max
        SignPass=countsp_max(tau3)/2
        PNumN(tau3,sequence)=sum(SPVector(tau3,1:SignPass,sequence))
        PNumI(tau3,sequence)=sum(SPVector(tau3,SignPass+1:countsp_max(tau3),sequence))
    end do
    end do

    do tau3=1, tau3_max
        serialnum=countsp_max(tau3)/2
        do sequence=1, TotalStates
            do count=1, serialnum
                ISPVector(tau3,count,sequence)=SPVector(tau3,count+serialnum,sequence)
                ISPVector(tau3,count+serialnum,sequence)=SPVector(tau3,count,sequence)
            end do
        end do
    end do

    ! do sequence=1, TotalStates
    !     write(*,*) sequence
    !     write(*,'(100I2)') SPVector(1,1:countsp_max(1),sequence), SPVector(2,1:countsp_max(2),sequence)
    !     write(*,'(100I2)') ISPVector(1,1:countsp_max(1),sequence), ISPVector(2,1:countsp_max(2),sequence)
    ! end do

    ! stop

    deallocate(ParticleVector)
    deallocate(IPhiVector)
    deallocate(PhiVector)
    deallocate(StateRemained)

end subroutine ParticleSort