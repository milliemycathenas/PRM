! subroutine ParticleSort

!     use ProcessVariables
!     use Constants
!     use PhysicalQuantities
!     use FunctionParameters

!     implicit none

!     integer :: OShellNum !number of single shells that contains odd number of particle
!     integer :: PNCStates !number of states with correct particle number on each j shell
!     integer :: EQStates, NEStates !number of states in EQVector or NEVector
!     integer, allocatable, dimension(:,:) :: PhiNIPair
!     integer, allocatable, dimension(:,:,:) :: NIPair
!     !serial number of a certain pair of particle states
!     !NIPair(tau3,1,sequence)=normal, NIPair(tau3,2,sequence)=inverse 
!     integer, allocatable, dimension(:,:) :: EQVector, NEVector
!     !EQVector stores states which are the same as their inverse states
!     !NEVector stores states which are different from their inverse states

!     if (.not. allocated(StateRemained)) allocate(StateRemained(DIM))
!     if (.not. allocated(PhiVector)) allocate(PhiVector(100,DIM))
!     if (.not. allocated(IPhiVector)) allocate(IPhiVector(100,DIM))
!     if (.not. allocated(ParticleVector)) allocate(ParticleVector(DIM))
!     if (.not. allocated(NIPair)) allocate(NIPair(10,2,DIM))
!     if (.not. allocated(PhiNIPair)) allocate(PhiNIPair(2,DIM))
!     if (.not. allocated(EQVector)) allocate(EQVector(sum(countsp_max),DIM))
!     if (.not. allocated(NEVector)) allocate(NEVector(sum(countsp_max),DIM))


!     OShellNum=0
!     do tau3=1, tau3_max
!         if (mod(Tau3PNum(tau3),2)/=0) OShellNum=OShellNum+1
!     end do

!     if (mod(OShellNum,2)/=0) then
        
!         do tau3=1, tau3_max
!             ParticleVector=0
!             StateRemained(tau3)=0

!             SortOccupationloop: &
!             do sequence=1, StateSequence(tau3)
!                 serialnum=sequence-1

!                 BinaryTransformloop: &
!                 do count=countsp_max(tau3), 1, -1
!                     ParticleVector(count) = mod(serialnum,2)
!                     serialnum=serialnum/2
!                 end do BinaryTransformloop

!                 if (sum(ParticleVector) == Tau3PNum(tau3)) then
!                     StateRemained(tau3)=StateRemained(tau3)+1

!                     do count=1, countsp_max(tau3)
!                         SPVector(tau3,count,StateRemained(tau3))=ParticleVector(count)
!                     end do
!                 end if
!             end do SortOccupationloop
!         end do

!         do tau3=1, tau3_max
!             serialnum=countsp_max(tau3)/2
!             do sequence=1, StateRemained(tau3)
!                 do count=1, serialnum
!                     ISPVector(tau3,count,sequence)=SPVector(tau3,count+serialnum,sequence)
!                     ISPVector(tau3,count+serialnum,sequence)=SPVector(tau3,count,sequence)
!                 end do
!             end do
!         end do

!         SignPass=1
!         do tau3=1, tau3_max
!             if (mod(Tau3PNum(tau3),2)/=0) then
!                 SignPass=tau3
!                 goto 100
!             end if
!         end do

!         100 do sequence=1, StateRemained(SignPass)
!                 NIPair(SignPass,1,sequence)=sequence
!             end do

!         NIPair(SignPass,2,:)=0
!         do column=1, StateRemained(SignPass)
!         do row=1, StateRemained(SignPass)
!             PassMark=1
!             do count=1, countsp_max(SignPass)
!                 if (ISPVector(SignPass,count,row)/=SPVector(SignPass,count,column)) PassMark=0
!             end do

!             if (PassMark==1) then
!                 NIPair(SignPass,2,column)=row
!             end if
!         end do
!         end do

!         serialnum=0
!         NStateCoupled=0
!         do sequence=1, StateRemained(SignPass)/2
!             if (NIPair(SignPass,1,sequence)/=0) then
!                 NStateCoupled(SignPass,sequence)=sequence
!                 NIPair(SignPass,1,sequence)=0
!                 NIPair(SignPass,1,NIPair(SignPass,2,sequence))=0
!                 serialnum=serialnum+1
!             end if
!         end do

!         StateNormal(SignPass)=serialnum

!         BinaryNumber=1
!         do tau3=1, tau3_max
!             BinaryNumber=BinaryNumber*StateSequence(tau3)
!         end do

!         PhiVector=0
!         TotalStates=0
!         do sequence=1, BinaryNumber

!             ParticleVector=0
!             serialnum=sequence-1
!             do count=sum(countsp_max), 1, -1
!                 ParticleVector(count) = mod(serialnum,2)
!                 serialnum=serialnum/2
!             end do

!             mark=1
!             do tau3=1, tau3_max
!                 EndPoint=sum(countsp_max(1:tau3))
!                 BeginPoint=EndPoint-countsp_max(tau3)+1

!                 if (sum(ParticleVector(BeginPoint:EndPoint))/=Tau3PNum(tau3)) mark=0
!             end do

!             if(mark==1) then

!                 index=0
!                 do label=1, StateNormal(SignPass)
!                     mark=1
!                     EndPoint=sum(countsp_max(1:SignPass))
!                     BeginPoint=EndPoint-countsp_max(SignPass)
!                     do count=1, countsp_max(SignPass)
!                         BeginPoint=BeginPoint+1
!                         if (ParticleVector(BeginPoint) /= &
!                         ISPVector(SignPass,count,NStateCoupled(SignPass,label))) mark=0
!                     end do
!                     if (mark==1) index=1
!                 end do

!                 if (index==0) then
                
!                     TotalStates=TotalStates+1

!                     do count=1, sum(countsp_max)
!                         PhiVector(count,TotalStates)=ParticleVector(count)
!                     end do
!                 end if
!             end if
!         end do

!     end if

!     if (mod(OShellNum,2)==0) then

!         if (OShellNum/=0) then

!             do tau3=1, tau3_max
!                 ParticleVector=0
!                 StateRemained(tau3)=0

!                 do sequence=1, StateSequence(tau3)
!                     serialnum=sequence-1

!                     do count=countsp_max(tau3), 1, -1
!                         ParticleVector(count) = mod(serialnum,2)
!                         serialnum=serialnum/2
!                     end do

!                     SignPass=countsp_max(tau3)/2
!                     PassMark=countsp_max(tau3)
!                     if (sum(ParticleVector(1:SignPass)) == Tau3PNum(tau3) .and. &
!                         sum(ParticleVector(SignPass+1:PassMark))==0) then
!                         StateRemained(tau3)=StateRemained(tau3)+1

!                         do count=1, countsp_max(tau3)
!                             SPVector(tau3,count,StateRemained(tau3))=ParticleVector(count)
!                         end do
!                     end if
!                 end do
!             end do

!             do tau3=1, tau3_max
!                 serialnum=countsp_max(tau3)/2
!                 do sequence=1, StateRemained(tau3)
!                     do count=1, serialnum
!                         ISPVector(tau3,count,sequence)=SPVector(tau3,count+serialnum,sequence)
!                         ISPVector(tau3,count+serialnum,sequence)=SPVector(tau3,count,sequence)
!                     end do
!                 end do
!             end do

!             TotalStates=0
!             do column=1, StateRemained(1)
!             do row=1, StateRemained(2)
!                 TotalStates=TotalStates+1
!                 JShellIndex(1,TotalStates)=column
!                 JShellIndex(2,TotalStates)=row
!             end do
!             end do

!             do sequence=1, TotalStates
!                 do tau3=1, tau3_max
!                     EndPoint=sum(countsp_max(1:tau3))
!                     BeginPoint=EndPoint-countsp_max(tau3)
!                 do count=1, countsp_max(tau3)
!                     PhiVector(count+BeginPoint,sequence)= &
!                     SPVector(tau3,count,JShellIndex(tau3,sequence))
!                 end do
!                 end do
!             end do

!             PassMark=TotalStates

!             do sequence=1, PassMark/2
!                 TotalStates=TotalStates+1

!                 EndPoint=sum(countsp_max(1:1))
!                 BeginPoint=EndPoint-countsp_max(1)
!                 do count=1, countsp_max(1)
!                     PhiVector(count+BeginPoint,TotalStates)= &
!                     SPVector(1,count,JShellIndex(1,sequence))
!                 end do
!                 EndPoint=sum(countsp_max(1:2))
!                 BeginPoint=EndPoint-countsp_max(2)
!                 do count=1, countsp_max(2)
!                     PhiVector(count+BeginPoint,TotalStates)= &
!                     ISPVector(2,count,JShellIndex(2,sequence))
!                 end do
!             end do

!             do sequence=1, PassMark/2
!                 TotalStates=TotalStates+1

!                 EndPoint=sum(countsp_max(1:1))
!                 BeginPoint=EndPoint-countsp_max(1)
!                 do count=1, countsp_max(1)
!                     PhiVector(count+BeginPoint,TotalStates)= &
!                     ISPVector(1,count,JShellIndex(1,sequence))
!                 end do
!                 EndPoint=sum(countsp_max(1:2))
!                 BeginPoint=EndPoint-countsp_max(2)
!                 do count=1, countsp_max(2)
!                     PhiVector(count+BeginPoint,TotalStates)= &
!                     SPVector(2,count,JShellIndex(2,sequence))
!                 end do
!             end do

!             ! do sequence=1, TotalStates
!             !     write(*,'(100I2)') PhiVector(1:sum(countsp_max),sequence)
!             ! end do

!         end if

!         if (OShellNum==0) then

!         end if 

!     end if

!     SPVector=0
!     ISPVector=0

!     do sequence=1, TotalStates
!         do tau3=1, tau3_max
!         EndPoint=sum(countsp_max(1:tau3-1))

!             do count=1, countsp_max(tau3)
!                 SPVector(tau3,count,sequence)=PhiVector(count+EndPoint,sequence)
!             end do
!         end do
!     end do

!     do sequence=1, TotalStates
!         do tau3=1, tau3_max
!             SignPass=countsp_max(tau3)/2
!             do count=1, countsp_max(tau3)
!                 ISPVector(tau3,count,sequence)=SPVector(tau3,count+SignPass,sequence)
!                 ISPVector(tau3,count+SignPass,sequence)=SPVector(tau3,count,sequence)
!             end do
!         end do
!     end do

!     do sequence=1, TotalStates
!     do tau3=1, tau3_max
!         SignPass=countsp_max(tau3)/2
!         PNumN(tau3,sequence)=sum(SPVector(tau3,1:SignPass,sequence))
!         PNumI(tau3,sequence)=sum(SPVector(tau3,SignPass+1:countsp_max(tau3),sequence))
!     end do
!     end do

!     deallocate(NEVector)
!     deallocate(EQVector)
!     deallocate(PhiNIPair)
!     deallocate(NIPair)
!     deallocate(ParticleVector)
!     deallocate(IPhiVector)
!     deallocate(PhiVector)
!     deallocate(StateRemained)

! end subroutine ParticleSort