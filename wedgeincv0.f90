!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   Objective: Write in a file the wavefunction scattered by a hard wall C
!   Compile  : gfortran wedgearb.f90 -O2 -o a.x -llapack -lblas          C
!   Files produced: geometria.dat -> contains the points of the barrier  C
!                   matriztkXX.Xx.dat -> contains the squared absolute   C
! values of the T-matrix elements                                        C
!                   psibwmkXX.XX.dat -> contains the denstity            C
! probability of the scattered wave function                             C
!                   varphikXX.XX.dat -> denstity                         C
! probability of the incident wave function                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

program  wedgeincv0

  REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)

    INTEGER, PARAMETER :: nr=438,n=300
    REAL(KIND=8), PARAMETER :: k = 100d0
    REAL(KIND=8), PARAMETER :: alfa = pi/12d0
    INTEGER, PARAMETER :: nb = 20
    REAL(KIND=8), PARAMETER :: R = DSIN(pi/3.0d0)/DSIN(2.0d0*pi/3.0d0-DSQRT(2.0d0))
    REAL(KIND=8), PARAMETER :: tetac = pi/3.0d0

    INTEGER i,j,aux1,aux2,imax
    REAL(KIND=8) :: ds
    REAL(KIND=8) :: xa(nr),ya(nr),xs(n,n),ys(n,n),lp,lcx
    COMPLEX(KIND=8) :: m(nr,nr),psi(n,n),varphi(n,n)
    CHARACTER(len=6) :: fileid1t
    CHARACTER(len=6) :: fileid2t

    aux1 = INT(k)
    aux2 = INT((k-aux1)*1000)
    WRITE(fileid1t,'(i0)') aux1
    WRITE(fileid2t,'(i0)') aux2

    OPEN(1,file='resultados/geometria.dat',STATUS='UNKNOWN')
    OPEN(2,file='resultados/matriztk'//trim(fileid1t)//'.'//trim(fileid2t)//'.dat',STATUS='UNKNOWN')
    OPEN(3,file='resultados/psibwmk'//trim(fileid1t)//'.'//trim(fileid2t)//'.dat',STATUS='UNKNOWN')
    OPEN(4,file='resultados/espaco.dat',STATUS='UNKNOWN')
    OPEN(8,file='resultados/varphi'//trim(fileid1t)//'.'//trim(fileid2t)//'.dat',STATUS='UNKNOWN')

    lp  = R
    ds  = lp/nr
    lcx = 2d0



    
    CALL DADOS(R,tetac,alfa,nr,n,nb,k)
    write(*,*) "Data ok!"

    CALL WALL(xa,ya,lp,ds,nr)
    DO i=1,nr
       WRITE(1,*) xa(i),ya(i)
    ENDDO
    CLOSE(1)
    write(*,*) "Wall ok!"

    
    Call MATRIX2(m,k,ds,xa,ya,nr,tetac)
    DO i=1,nr
       DO j=1,nr
          WRITE(2,*) i,j,REAL(m(i,j)*DCONJG((m(i,j))))
       ENDDO
       WRITE(2,*)"                                           "
    ENDDO
    CLOSE(2)
    write(*,*) "Matrix ok!"

    CALL SPACE(xs,ys,lcx,tetac,n)
    DO i=1,n
       DO j=1,n
          WRITE(4,*) xs(i,j),ys(i,j),SQRT(i*j*1.0d0)
       ENDDO
       WRITE(4,*)"                                    "
    ENDDO
    CLOSE(4)
    write(*,*) "Space ok!"

    
    CALL WAVEFUNCTION2(k,m,xa,ya,xs,ys,ds,psi,nr,n,nb,tetac,alfa,varphi)
    DO i=1,n
      imax = min(n,INT(DTAN(tetac)*i))
       DO j=1,n
          if ( j.LE.imax ) then
            WRITE(3,*) xs(i,j),ys(i,j),REAL(psi(i,j)*DCONJG((psi(i,j))))
            WRITE(8,*) xs(i,j),ys(i,j), REAL(varphi(i,j)*DCONJG((varphi(i,j))))
          else
            WRITE(3,*) xs(i,j),ys(i,j),0d0
            WRITE(8,*) xs(i,j),ys(i,j),0d0
          end if
       ENDDO
       WRITE(3,*)"                                           "
       WRITE(8,*)"                                           "
    ENDDO
    CLOSE(3)
    close(8)
    write(*,*) "wavefunction ok!"
    STOP
  END program  wedgeincv0
  
  SUBROUTINE WALL(x,y,lp,ds,nr)
    REAL(KIND=8) :: x(nr),y(nr),lp,ds
    REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)
    INTEGER :: nr,i


    DO i=1,nr

      x(i) = 1d0-DCOS(DSQRT(2d0))*(DBLE(i)*ds-ds/2.0d0)
      y(i) = DSIN(DSQRT(2d0))*(DBLE(i)*ds-ds/2.0d0)
    ENDDO

    RETURN
  END SUBROUTINE WALL



  SUBROUTINE SPACE(xs,ys,lcx,tetac,n)
    REAL(KIND=8) :: xs(n,n),ys(n,n),lcx,lcy,tetac,dsx,dsy
    INTEGER      :: n,i,imax
    lcy = lcx
    dsx = lcx/n
    dsy = lcy/n
    DO i=1,n
       DO j=1,n
        xs(i,j) = DBLE(i)*dsx-dsx/2.0d0
        ys(i,j) = DBLE(j)*dsy-dsy/2.0d0
       ENDDO
    ENDDO
  END SUBROUTINE SPACE


  SUBROUTINE dgreta(ds,k,diag)
    COMPLEX(KIND=8) :: diag,auxc,int
    REAL(KIND=8) :: den,parc,y1,y2,pi,ds,k
    INTEGER :: n
    REAL(KIND=8), PARAMETER :: euler=0.5772156649015329d0
    pi = DACOS(-1.d0)
    n = 0
    den = 1.d0
    parc = 0.d0
    int = (0.d0,0.d0)
  700 CONTINUE
    den = den*DBLE(2*n)**2
    IF(n==0) den = 1.d0
    y1 = (2.d0/(2.d0*DBLE(n)+1.d0))*(ds/2.d0)**(2*n+1)*k**(2*n)/den
    IF(n/=0) parc = parc + 1.d0/DBLE(n)
    y2 = (2.d0/pi)*y1*(euler-parc-1.d0/(2.d0*DBLE(n)+1.d0)+DLOG(k*ds/4.d0))
    auxc = DCMPLX(y1,y2)
    int = int + (-1.d0)**n*(auxc)
    n = n + 1
    IF(n>=100) THEN
       WRITE(*,*) 'Muitas iteracoes para a diagonal 1'
       WRITE(*,*) y1,y2,int
       STOP
    ENDIF
    IF(DABS(y1)>1.d-20.OR.DABS(y2)>1.d-20) GOTO 700
    diag = int
    RETURN
  END SUBROUTINE dgreta

  SUBROUTINE DADOS(lado,angcunha,alfa,nr,n,nb,k)
    INTEGER nr,n,aux1,aux2,nb
    REAL(KIND=8) :: lado,angcunha,alfa,k
    CHARACTER(len=6) :: fileid4t
    CHARACTER(len=6) :: fileid5t

    aux1 = INT(k)
    aux2 = INT((k-aux1)*1000)
    WRITE(fileid4t,'(i0)') aux1
    WRITE(fileid5t,'(i0)') aux2
    OPEN(unit=10,file='resultados/dados'//trim(fileid4t)//'.'//trim(fileid5t)//'.dat',STATUS='UNKNOWN')
    WRITE(10,*) 'LAdo', lado
    WRITE(10,*) 'angcunha', angcunha
    WRITE(10,*) 'alfa', alfa
    WRITE(10,*) 'k', k
    WRITE(10,*) 'nr', nr
    WRITE(10,*) 'n', n
    WRITE(10,*) 'nb', nb
    CLOSE(10)

  END SUBROUTINE DADOS


  SUBROUTINE MATRIX2(m,k,ds,xa,ya,nr,angcunha)
      COMPLEX(KIND=8), PARAMETER :: sigma=(0d0,-0.25d0)
      REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)
      INTEGER :: i,j,nr,ipiv(nr),info,l,Maiusculo
      REAL(KIND=8) :: k,xa(nr),ya(nr),y1,y2,ds,x1
      REAL(KIND=8) :: raux(2),angcunha
      COMPLEX(KIND=8) :: intg,diag,G0,m(nr,nr),work(nr,nr),G1

      CALL dgreta(ds,k,diag)

      Maiusculo = Int(pi/angcunha)
      DO i=1,nr
         DO j=1,i

            IF(i==j)THEN
              G1 = diag/ds
              G0 = (0.0d0,0.0d0)
              do l = 1, 2*Maiusculo-1
                raux=(0.0d0,0.0d0)
                Call imagem(xa(j),ya(j),l,angcunha,raux)
                x1 = k*DSQRT((raux(1)-xa(i))**2+(raux(2)-ya(i))**2)
                y1   = bessel_j0(x1)
                y2   = bessel_y0(x1)
                G0   = G0 + ((-1.0d0)**l)*DCMPLX(y1,y2)
              end do
              G0= G0 + G1

            ELSE
              G0 = (0.0d0,0.0d0)
              do l = 0, 2*Maiusculo-1
                raux=(0.0d0,0.0d0)
                Call imagem(xa(j),ya(j),l,angcunha,raux)

                x1 = k*DSQRT((raux(1)-xa(i))**2+(raux(2)-ya(i))**2)
                y1   = bessel_j0(x1)
                y2   = bessel_y0(x1)
                G0   = G0 + ((-1.0d0)**l)*DCMPLX(y1,y2)

              end do

            ENDIF
            intg   = sigma*G0*ds
            m(i,j) = intg
            m(j,i) = intg
            IF(CDABS(m(i,j))==0.d0)THEN
               WRITE(*,*) i,j,m(i,j),x1
               WRITE(*,*) intg, ds, G0
               stop 'Problemas'
            ENDIF
         ENDDO
      ENDDO

      CALL ZGETRF(nr,nr,m,nr,ipiv,info)
      CALL ZGETRI(nr,m,nr,ipiv,work,nr,info)
      RETURN
    END SUBROUTINE matrix2


    SUBROUTINE imagem(x,y,l,angcunha,result)
    REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)
    REAL(KIND=8) :: x,y,angcunha, ref(2,2),rot(2,2), aux(2),result(2)
    INTEGER      :: l

    aux(1) = x
    aux(2) = y

    if ( Mod(l,2)==0 ) then
      rot(1,1)= DCos(angcunha*DBLE(l))
      rot(2,2)= rot(1,1)
      rot(1,2)= -DSin(angcunha*DBLE(l))
      rot(2,1)= -rot(1,2)
      result = matmul(rot,aux)
    else
      ref(1,1)= DCos(angcunha*DBLE(l+1))
      ref(2,2)= -ref(1,1)
      ref(1,2)= DSin(angcunha*DBLE(l+1))
      ref(2,1)= ref(1,2)
      result = matmul(ref,aux)
    end if
    RETURN
  END SUBROUTINE imagem

  SUBROUTINE WAVEFUNCTION2(k,m,xa,ya,xs,ys,ds,psi,nr,n,nb,angcunha,alfa,varphi)

      COMPLEX(KIND=8), PARAMETER :: sigma=(0d0,-0.25d0)
      REAL(KIND=8), PARAMETER :: pi=DACOS(-1.0d0)
      INTEGER :: i,j,l,nr,n,nb,ii,Maiusculo,imax
      REAL(KIND=8) :: k,xs(n,n),ys(n,n),xa(nr),ya(nr),y1,y2,angcunha,alfa
      REAL(KIND=8) :: ds,x1,raux(2)
      COMPLEX(KIND=8) :: intg,G0,m(nr,nr),psi(n,n),tphi(nr),phi,varphi(n,n)
      COMPLEX(KIND=8), PARAMETER :: ui=(0.0d0,1.0d0)

      Maiusculo = Int(pi/angcunha)
      DO i=1,nr
         tphi(i)=(0.d0,0.d0)
         DO j=1,nr
            phi=0
            do ll = 0, 2*Maiusculo-1
              raux=(0.0d0,0.0d0)
            
              Call imagem(xa(j),ya(j),ll,angcunha,raux)
              
              x1 = k*(raux(1)*Dcos(alfa+pi)+raux(2)*DSIN(alfa+pi))
              phi   = phi + ((-1.0d0)**ll)*((1d0,0)*DCos(x1)+ui*Dsin(x1))
            end do
            tphi(i) = tphi(i)+m(i,j)*phi
         ENDDO
      ENDDO

      DO i=1,n
        imax =min(n,INT(DTAN(angcunha)*i))
          DO j=1,imax
            intg = (0.d0,0.0d0)
            DO l=1,nr
              G0 = (0.0d0,0.0d0)
              do ll = 0, 2*Maiusculo-1
                raux=(0.0d0,0.0d0)
                Call imagem(xa(l),ya(l),ll,angcunha,raux)
                x1 = k*DSQRT((raux(1)-xs(i,j))**2+(raux(2)-ys(i,j))**2)
                y1   = bessel_j0(x1)
                y2   = bessel_y0(x1)
                G0   = G0 + ((-1.0d0)**ll)*DCMPLX(y1,y2)
              end do
               G0   = sigma*G0
               intg = intg+G0*tphi(l)*ds
            ENDDO


            phi=0
            do ll = 0, 2*Maiusculo-1
              raux=(0.0d0,0.0d0)
            
              Call imagem(xs(i,j),ys(i,j),ll,angcunha,raux)
              
              x1 = k*(raux(1)*Dcos(alfa+pi)+raux(2)*DSIN(alfa+pi))
              phi   = phi + ((-1.0d0)**ll)*((1d0,0)*DCos(x1)+ui*Dsin(x1))
            end do
            varphi(i,j) = phi
            psi(i,j) = phi-intg
          ENDDO
          Write(*,*) "i = ",i
      ENDDO
      RETURN
    END SUBROUTINE wavefunction2