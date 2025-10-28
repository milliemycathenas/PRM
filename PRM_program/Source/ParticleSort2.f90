! subroutine ParticleSort

!     use ProcessVariables
!     use Constants
!     use PhysicalQuantities
!     use FunctionParameters

!     implicit none

!     if (.not. allocated(StateRemained)) allocate(StateRemained(DIM))
!     if (.not. allocated(PhiVector)) allocate(PhiVector(100,DIM))
!     if (.not. allocated(IPhiVector)) allocate(IPhiVector(100,DIM))
!     if (.not. allocated(ParticleVector)) allocate(ParticleVector(DIM))
!     if (.not. allocated(NPassArray)) allocate(NPassArray(10,DIM))
!     if (.not. allocated(IPassArray)) allocate(IPassArray(10,DIM))

!     do tau3=1, tau3_max

!         ParticleVector=0
!         StateRemained(tau3)=0

!         SortOccupationloop: &
!         do sequence=1, StateSequence(tau3)
!             serialnum=sequence-1

!             BinaryTransformloop: &
!             do count=countsp_max(tau3), 1, -1
!                 ParticleVector(count) = mod(serialnum,2)
!                 serialnum=serialnum/2
!             end do BinaryTransformloop

!             if (sum(ParticleVector) == Tau3PNum(tau3)) then
!                 StateRemained(tau3)=StateRemained(tau3)+1

!                 do count=1, countsp_max(tau3)
!                     SPVector(tau3,count,StateRemained(tau3))=ParticleVector(count)
!                 end do
!             end if
!         end do SortOccupationloop
!     end do

!     do tau3=1, tau3_max
!         serialnum=countsp_max(tau3)/2
!         do sequence=1, StateRemained(tau3)
!             do count=1, serialnum
!                 ISPVector(tau3,count+serialnum,sequence)=SPVector(tau3,count,sequence)
!             end do
!             do count=1, serialnum
!                 ISPVector(tau3,count,sequence)=SPVector(tau3,count+serialnum,sequence)
!             end do
!         end do
!     end do

!     do tau3=1, tau3_max
!         do sequence=1, StateRemained(tau3)
!             NPassArray(tau3,sequence)=sequence
!             IPassArray(tau3,sequence)=sequence
!         end do
!     end do

!     do tau3=1, tau3_max
!         sequence=0
!         do column=1, StateRemained(tau3)
!         do row=1, StateRemained(tau3)
!             SignPass=1
!             do count=1, countsp_max(tau3)
!                 if (ISPVector(tau3,count,row)/=SPVector(tau3,count,column)) SignPass=0
!             end do

!             if (SignPass==1) then
!                 if (IPassArray(tau3,column)/=0) then
!                     sequence=sequence+1
!                     IPassArray(tau3,row)=0
!                     NStateCoupled(tau3,sequence)=row
!                 end if
!             end if
!         end do
!         end do
!         StateNormal(tau3)=sequence
!     end do

!     do tau3=1, tau3_max
!         sequence=0
!         do serialnum=1, StateRemained(tau3)
!             if (IPassArray(tau3,serialnum)/=0) then
!                 sequence=sequence+1
!                 IStateCoupled(tau3,sequence)=IPassArray(tau3,serialnum)
!             end if
!         end do
!     end do

!     BinaryNumber=1
!     do tau3=1, tau3_max
!         BinaryNumber=BinaryNumber*StateSequence(tau3)
!     end do

!     SignPass=1
!     do tau3=1, tau3_max
!         if (mod(Tau3PNum(tau3),2)/=0) then
!             SignPass=tau3
!             ! goto 100
!         end if
!     end do
    
!     100 PhiVector=0
!     TotalStates=0
!     do sequence=1, BinaryNumber

!         ParticleVector=0
!         serialnum=sequence-1
!         do count=sum(countsp_max), 1, -1
!             ParticleVector(count) = mod(serialnum,2)
!             serialnum=serialnum/2
!         end do

!         mark=1
!         do tau3=1, tau3_max
!             EndPoint=sum(countsp_max(1:tau3))
!             BeginPoint=EndPoint-countsp_max(tau3)+1

!             if (sum(ParticleVector(BeginPoint:EndPoint))/=Tau3PNum(tau3)) mark=0
!         end do

!         if(mark==1) then

!             index=0
!             do label=1, StateNormal(SignPass)
!                 mark=1
!                 EndPoint=sum(countsp_max(1:SignPass))
!                 BeginPoint=EndPoint-countsp_max(SignPass)
!                 do count=1, countsp_max(SignPass)
!                     BeginPoint=BeginPoint+1
!                     if (ParticleVector(BeginPoint) /= &
!                     ISPVector(SignPass,count,NStateCoupled(SignPass,label))) mark=0
!                 end do
!                 if (mark==1) index=1
!             end do

!             if (index==0) then
            
!                 TotalStates=TotalStates+1

!                 do count=1, sum(countsp_max)
!                     PhiVector(count,TotalStates)=ParticleVector(count)
!                 end do
!             end if
!         end if
!     end do

!     SPVector=0
!     ISPVector=0

!     do sequence=1, TotalStates
!         do tau3=1, tau3_max
!             EndPoint=sum(countsp_max(1:tau3))
!             BeginPoint=EndPoint-countsp_max(tau3)
!             do count=1, countsp_max(tau3)
!                 BeginPoint=BeginPoint+1
!                 SPVector(tau3,count,sequence)=PhiVector(BeginPoint,sequence)
!             end do
!             PNumN(tau3,sequence)=sum(SPVector(tau3,1:countsp_max(tau3)/2,sequence))
!             PNumI(tau3,sequence)=sum(SPVector(tau3,countsp_max(tau3)/2+1:countsp_max(tau3),sequence))
!         end do
!     end do

!     do tau3=1, tau3_max
!         serialnum=countsp_max(tau3)/2
!         do sequence=1, TotalStates
!             do count=1, serialnum
!                 ISPVector(tau3,count+serialnum,sequence)=SPVector(tau3,count,sequence)
!             end do
!             do count=1, serialnum
!                 ISPVector(tau3,count,sequence)=SPVector(tau3,count+serialnum,sequence)
!             end do
!         end do
!     end do

!     deallocate(NPassArray)
!     deallocate(IPassArray)
!     deallocate(ParticleVector)
!     deallocate(IPhiVector)
!     deallocate(PhiVector)
!     deallocate(StateRemained)

! end subroutine ParticleSort