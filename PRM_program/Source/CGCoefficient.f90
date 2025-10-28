subroutine CGCoefficient(j1, j2, j3, m1, m2, m3)

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use CGParameters

    implicit none 

    real :: FactorialRes

    real :: j1, j2, j3 !angular momentum of bra&ket, j3 is the total angular momentum
    real :: m1, m2, m3 !magnetic quantum number of bra&ket, m3 is the total magnetic quantum number
    !real :: CGResult

    CGResult = 0

    SignPass = 1
    if (m1+m2/=m3) SignPass = 0
    if (j1+j2-j3<0) SignPass = 0
    if (j2+j3-j1<0) SignPass = 0
    if (j3+j1-j2<0) SignPass = 0
    if (j1+j2+j3+1<0) SignPass = 0
    if (j1+m1<0) SignPass = 0
    if (j1-m1<0) SignPass = 0
    if (j2+m2<0) SignPass = 0
    if (j2-m2<0) SignPass = 0
    if (j3+m3<0) SignPass = 0
    if (j3-m3<0) SignPass = 0

    if ( SignPass==1 ) then

        CGResult = 2*j3+1
        CGResult = CGResult * FactorialRes(int(j1+j2-j3))*FactorialRes(int(j2+j3-j1))*FactorialRes(int(j3+j1-j2)) &
        *FactorialRes(int(j1-m1))*FactorialRes(int(j1+m1)) &
        *FactorialRes(int(j2-m2))*FactorialRes(int(j2+m2)) &
        *FactorialRes(int(j3-m3))*FactorialRes(int(j3+m3)) &
        /FactorialRes(int(j1+j2+j3+1))
        CGResult = sqrt(CGResult)

        SumMuRes = 0
        do serialnum=0, 10

            Factorial123 = int(j1+j2-j3-serialnum)
            Factorial11 = int(j1-m1-serialnum)
            Factorial22 = int(j2+m2-serialnum)
            Factorial321 = int(j3-j2+m1+serialnum)
            Factorial312 = int(j3-j1-m2+serialnum)

            if (Factorial123>=0 .and. Factorial11>=0 .and. Factorial22>=0 &
            .and. Factorial321>=0 .and. Factorial312>=0) then 
                ! write(*,*) FactorialRes(serialnum)* FactorialRes(Factorial123)* FactorialRes(Factorial11)* &
                ! FactorialRes(Factorial22)*FactorialRes(Factorial321)* FactorialRes(Factorial312)

                SumMuRes = SumMuRes + (1/real(((-1)**(serialnum))*(FactorialRes(serialnum) &
                *FactorialRes(Factorial123)*FactorialRes(Factorial11)*FactorialRes(Factorial22)&
                *FactorialRes(Factorial321)*FactorialRes(Factorial312))))
            end if

        end do

        CGResult = CGResult*SumMuRes
    
    else 
        CGResult = 0
    end if


end subroutine CGCoefficient