subroutine ReAlign

	use ProcessVariables
	use Constants
	use PhysicalQuantities
	use DiagonalVariables
	use FunctionParameters

	implicit none

    do tau3=1, tau3_max
        PassMark=countsp_max(tau3)/2
        do count=1, PassMark
            EAlign(tau3,2,count)=EAlign(tau3,1,count)
            EAlign(tau3,2,count+PassMark)=EAlign(tau3,1,count)
        end do
    end do 

    do tau3=1, tau3_max
        PassMark=countsp_max(tau3)/2
        if ( mod(jsp(tau3)+0.5, 2.)==0 ) then
            do row=1, countsp_max(tau3)
            do column=1, PassMark
                C_nu(tau3,2*row-1,column)=CAlign(tau3,1,row,column)
                C_nu(tau3,2*row,column)=0
            end do
            do column=PassMark+1, countsp_max(tau3)
                C_nu(tau3,2*row-1,column)=0
                C_nu(tau3,2*row,column)=((-1)**(NilVector(tau3,1,row)-NilVector(tau3,2,row))) &
                *CAlign(tau3,1,PassMark+1-row,column-PassMark)
            end do
            end do
        else
            do row=1, countsp_max(tau3)
            do column=1,PassMark
                C_nu(tau3,2*row,column)=CAlign(tau3,1,row,column)
                C_nu(tau3,2*row-1,column)=0
            end do
            do column=PassMark+1, countsp_max(tau3)
                C_nu(tau3,2*row,column)=0
                C_nu(tau3,2*row-1,column)=((-1)**(NilVector(tau3,1,row)-NilVector(tau3,2,row))) &
                *CAlign(tau3,1,PassMark+1-row,column-PassMark)
            end do
            end do
        end if
    end do

    ! do tau3=1, tau3_max
    !     write(*,*) 'tau3', tau3
    !     do row=1, countsp_max(tau3)
    !         write(*,'(100f12.6)') C_nu(tau3,row,1:countsp_max(tau3))
    !     end do
    ! end do

    ! stop
    
end subroutine ReAlign