subroutine VInitialise

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters
    use CGParameters

    implicit none

    NJ3Mono = 0
    NJ3SMono = 0
    NJ3SPMono = 0
    NJ3STB = 0
    NJSMono = 0
    NJSPMono = 0
    NJSTB = 0
    NJuJdMono = 0
    NJuJdPMono = 0
    NJuJdTB = 0
    NJdJuMono = 0
    NJdJuPMono = 0
    NJdJuTB = 0
    NJuMono = 0
    NJdMono = 0

    IJ3Mono = 0   
    IJ3SMono = 0
    IJ3SPMono = 0
    IJ3STB = 0
    IJSMono = 0
    IJSPMono = 0
    IJSTB = 0
    IJuJdMono = 0
    IJuJdPMono = 0
    IJuJdTB = 0
    IJdJuMono = 0
    IJdJuPMono = 0
    IJdJuTB = 0
    IJuMono = 0
    IJdMono = 0

end subroutine VInitialise

subroutine RInitialise

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters

    implicit none

    RNJUDSMono=0
    RNJUDSPMono=0
    RNJUDSTB=0
    RNJxySMono=0
    RNJxySPMono=0
    RNJxySTB=0
    RNJUp=0
    RNJDown=0

    RNJ3=0
    RNJ3SMono=0
    RNJ3SPMono=0
    RNJ3STB=0

    RIJUDSMono=0
    RIJUDSPMono=0
    RIJUDSTB=0
    RIJxySMono=0
    RIJxySPMono=0
    RIJxySTB=0
    RIJUp=0
    RIJDown=0

    RIJ3=0
    RIJ3SMono=0
    RIJ3SPMono=0
    RIJ3STB=0

end subroutine RInitialise