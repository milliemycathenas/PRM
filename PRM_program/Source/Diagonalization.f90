    SUBROUTINE Diagonalization(N1,N,EV,A,D,V,IROT,B,Z)

                   !N1：所要求解的若干哈密顿量矩阵中的最大维数
                   !N：当前求解的哈密顿量矩阵的维数
                   !EV：逻辑量，取.TRUE.表示输出本征态矢，取.FALSE.表示不输出本征态矢
                   !该逻辑量在调用JACOBI的程序中也要加以说明;LOGICAL EV,可以换成一个整型量分别取1,0代替
                   !A：二维数组，输入哈密顿量矩阵元
                   !D：一维数组，输出本征值
                   !V：二维数组，输出本征态矢
                   !IROT：子程序中的一个输出参量，由具体算法而来，表示求解过程中哈密顿量变换的次数，可不必管它。
                   !B, Z：都是二维工作数组，运算过程中要用到，在调用它们的程序（可以是主程序，也可以说是某个子程序）中也要加以说明                   
    INTEGER EV
	!DIMENSION A(N1,N1),D(N1),V(N1,N1),B(N1),Z(N1)
    DIMENSION A(N1,N1)
    DIMENSION D(N1)
    DIMENSION V(N1,N1)
    DIMENSION B(N1)
    DIMENSION Z(N1)

    !DO I=1,N
    !    DO J=1,N
    !WRITE (*,*) 'A(',J,',',I,')',A(J,I)
    !    END DO
    !END DO
    !WRITE (*,*) '子程序里的求解维数',N,N1

	IF(EV==0) GOTO 10
	DO IP=1,N
	  DO IQ=1,N
	    IF(IP-IQ) 50,60,50
60	    V(IP,IP)=1.
	    GOTO 70
50	    V(IP,IQ)=0.
70      END DO
	END DO
10	DO IP=1,N
	  D(IP)=A(IP,IP)
	  B(IP)=D(IP)
	  Z(IP)=0.
	END DO
	IROT=0
	DO I=1,50
	  SM=0.
	  NM1=N-1
	  DO IP=1,NM1
	    IPP1=IP+1
	    DO IQ=IPP1,N
	      SM=SM+ABS(A(IP,IQ))
	    END DO
	  END DO
	  IF(SM) 110,120,110
110	  IF(I-4) 130,140,140
130	  TRESH=0.2*SM/(FLOAT(N)*FLOAT(N))
	  GOTO 150
140	  TRESH=0.
150	  DO IP=1,NM1
	    IPP1=IP+1
	    DO IQ=IPP1,N
	      G=100.*ABS(A(IP,IQ))
	      IF(I.GT.4.AND.ABS(D(IP))+G.EQ.ABS(D(IP))&
          .AND.ABS(D(IQ))+G.EQ.ABS(D(IQ))) GOTO 200
		  IF(ABS(A(IP,IQ)).LE.TRESH) GOTO 160
	      H=D(IQ)-D(IP)
	      IF(ABS(H)+G.EQ.ABS(H)) GOTO 240
	      THETA=0.5*H/A(IP,IQ)
	      T=1./(ABS(THETA)+SQRT(1.+THETA*THETA))
	      IF(THETA.LT.0.) T=-T
	      GOTO 250
240	      T=A(IP,IQ)/H
250	      C=1./SQRT(1.+T*T)
	      S=T*C
	      H=T*A(IP,IQ)
	      Z(IP)=Z(IP)-H
	      Z(IQ)=Z(IQ)+H
	      D(IP)=D(IP)-H
	      D(IQ)=D(IQ)+H
	      A(IP,IQ)=0.
	      IPM1=IP-1
	      IF(IPM1) 260,260,270
270	      DO J=1,IPM1
	        G=A(J,IP)
	        H=A(J,IQ)
	        A(J,IP)=C*G-S*H
	        A(J,IQ)=S*G+C*H
	      END DO
260	      IQM1=IQ-1
	      IF(IQM1-IPP1) 300,290,290
290	      DO J=IPP1,IQM1
	        G=A(IP,J)
	        H=A(J,IQ)
	        A(IP,J)=C*G-S*H
              A(J,IQ)=S*G+C*H
            END DO
300	      IQP1=IQ+1
	      IF(N-IQP1) 330,320,320
320	      DO J=IQP1,N
	        G=A(IP,J)
	        H=A(IQ,J)
	        A(IP,J)=C*G-S*H
	        A(IQ,J)=S*G+C*H
	      END DO
330	      IF(EV==0) GOTO 350
	      DO J=1,N
	        G=V(J,IP)
	        H=V(J,IQ)
	        V(J,IP)=C*G-S*H
	        V(J,IQ)=S*G+C*H
            END DO
350	      IROT=IROT+1
	      GOTO 160
200	      A(IP,IQ)=0.
160	    END DO
	  END DO
	  DO IP=1,N
	    B(IP)=B(IP)+Z(IP)
	    D(IP)=B(IP)
	    Z(IP)=0.
	  END DO
	END DO
120	RETURN
    END SUBROUTINE Diagonalization