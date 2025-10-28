    SUBROUTINE EnergyAlignment(N,M,NEV,E,C)
    !N: the max dimension
    !M: the dimension of the present matrix
    !NEV: whether output the engin state (1 or 0)
    !E: 1-dimension matrix for eigen-energy
    !C: 2-dimension matrix for process, only need to be defined
    
        INTEGER H
        LOGICAL INTCH
        DIMENSION E(N),C(N,N)
        H=M
1       IF(H.GT.1) THEN
          H=(H+1)/2
2         INTCH=.FALSE.
          MM=M-H
          DO K=1,MM
            L=K+H
            IF(E(K).GT.E(L)) THEN
	          TEMP=E(K)
	          E(K)=E(L)
	          E(L)=TEMP
	          IF(NEV.EQ.1) THEN
	            DO I=1,M
		          TEMP=C(I,K)
		          C(I,K)=C(I,L)
		          C(I,L)=TEMP
		        END DO
	          END IF
	          INTCH=.TRUE.
	        END IF
	      END DO
	      IF(INTCH) GO TO 2
	      GO TO 1
	    END IF
	END SUBROUTINE EnergyAlignment