subroutine CalDimension

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use CGParameters

    implicit none

    if (.not. allocated(Combination)) allocate(Combination(10))

    do tau3=1, tau3_max
        call CalCombination(countsp_max(tau3),Tau3PNum(tau3))
        Combination(tau3) = CombinationRes
    end do

    DimensionNumber = 1
    do tau3=1, tau3_max
        DimensionNumber = DimensionNumber*Combination(tau3)
    end do
    DimensionNumber = 0.25 * (2*I+1) * DimensionNumber

    deallocate(Combination)
    
end subroutine CalDimension

function FactorialRes(n)

    use ProcessVariables
    use CGParameters

    implicit none

    integer :: FactorialRes
    integer :: n !this factorial starts from n

    FactorialRes = n

    if (n<0) then
        write(*,*) "error: the factorial factor must be non-negative."
    else
        if ((n==0) .or. (n==1)) then
            FactorialRes = 1
        else
            do count=n-1, 1, -1
                FactorialRes = FactorialRes * count
            end do
        end if
    end if

    return

end function FactorialRes

subroutine CalCombination(n, m)

    use CGParameters

    implicit none

    integer :: FactorialRes
    integer :: n, m !for C_n^m
    integer :: factor_nm, factor_n, factor_m

    factor_nm = FactorialRes(n-m)

    factor_n = FactorialRes(n)
 
    factor_m = FactorialRes(m)

    CombinationRes = factor_n / ( factor_nm * factor_m )

end subroutine CalCombination