subroutine KPlus

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    if (NSCVector(1,bra)==NSCVector(1,ket)+1) then

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count) = SPVector(tau3,count,NSCVector(2,bra))
                KSPVector(tau3,count) = SPVector(tau3,count,NSCVector(2,ket))
            end do
        end do

        NPlusMloop: &
        do tau3=1, tau3_max
            NMMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle NMMuloop
            NMNuloop: &
            do AnnihNu=1, countsp_max(tau3) 
                if (KSPVector(tau3,AnnihNu)==0) cycle NMNuloop
                
                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence = sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence = sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    NJuMono = NJuMono+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    NJdMono = NJdMono+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do NMNuloop
            end do NMMuloop
        end do NPlusMloop

    end if

    if (NSCVector(1,bra)==ISCVector(1,ket)+1) then

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count) = SPVector(tau3,count,NSCVector(2,bra))
                KSPVector(tau3,count) = ISPVector(tau3,count,NSCVector(2,ket))
            end do
        end do

        IPlusMloop: &
        do tau3=1, tau3_max
            IMMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle IMMuloop
            IMNuloop: &
            do AnnihNu=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNu)==0) cycle IMNuloop

                BSPVector(tau3,AnnihMu) = BSPVector(tau3,AnnihMu)-1
                MuSequence = sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu) = KSPVector(tau3,AnnihNu)-1
                NuSequence = sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    IJuMono = IJuMono+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    IJdMono = IJdMono+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do IMNuloop
            end do IMMuloop
        end do IPlusMloop

    end if

end subroutine KPlus