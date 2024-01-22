!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   Objective: Write in a file the mean value of the squared module of   C
! the T-matrix elements for each k. This values are used to find the     C
! numeric resonant k close to the initial guess "kanalitico".            C
! We vary the search range of k to meet the desired precision            C
!   Compile  : gfortran wedgearbresv6.f90 -O2 -o a.x -llapack -lblas     C

program wedgearbresv6

  REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)
  REAL(KIND=8), PARAMETER :: pre = 0.01d0
  !razao
  REAL(KIND=8), PARAMETER :: r = 0.05d0
  REAL(KIND=8), PARAMETER :: tetac = pi/3.0d0 
  REAL(KIND=8), PARAMETER :: Kanalitico = 101.98d0

    INTEGER j,jkfim,ip,nr,nrmax
    REAL(KIND=8) :: ds,k,kini,kfim,dk,lp,lcx,kaux,kaux2,mediaaux,lado
    REAL(KIND=8) :: media, desviopadrao, errodamedia,max, kmax,dpmax,emmax
    COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE :: M
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: xa,ya

    OPEN(2,file='resultados/ress.dat',STATUS='UNKNOWN')
    OPEN(3,file='resultados/ressprec.dat',STATUS='UNKNOWN')
    OPEN(4,file='resultados/kressonante.dat',STATUS='UNKNOWN')

    lado = DSQRT(3d0)
    lp  = lado
    lcx = 3.0d0

    dk = 0.05d0
  kini =  Kanalitico - dk
  kfim = kini+2*dk

  kaux = 0
  kaux2 = 0
  mediaaux = 0
  do ip = 1, 200
    
    dk = (kfim-kini)/10 

    jkfim = INT((kfim-kini)/dk)
    kmax = 0.0d0
    max = 0.0d0
    do j = 1, jkfim
      k= kini + dk*(DBLE(j)-1)

      nr = INT(lp*k/(2*pi*r))

      ds  = lp/nr
      ALLOCATE(xa(1:nr),ya(1:nr),M(1:nr,1:nr))
      CALL WALL(xa,ya,lp,ds,nr)

      m=0.0d0
      Call MATRIX2(m,k,ds,xa,ya,nr,tetac)
      media = 0.0d0
      Call submedia(m,nr,media)
      desviopadrao = 0.0d0
      Call subdp(m,nr,media,desviopadrao)
      errodamedia= desviopadrao/Dble(nr)
      WRITE(2,*) k, media, desviopadrao, errodamedia, nr

      if ( MOD(j,5)==0 ) then
        write(*,*) "j=",j  
      end if 

      if ( media>max ) then
        max = media
        kmax = k
        dpmax = desviopadrao
        emmax = errodamedia
        nrmax = nr 
      end if
      DEALLOCATE(xa,ya,M)
    end do

    WRITE(3,*) kmax, max, dpmax, emmax, nrmax
    if ( max>mediaaux ) then
      kaux2 = kaux
      kaux = kmax
      mediaaux = max
    end if
    if ( DABS(kaux-kaux2)<pre) then
      Close(2)
      Close(3)
      write(*,*) "k ressonante =", kaux
      write(*,*) "nr           =", nr
      write(*,*) "DeltaK% =", DABS(kaux-Kanalitico)*100d0/Kanalitico
      write(4,*) "k ressonante =", kaux
      Close(4)
      STOP
    end if

    kini = kmax - 3*dk
    kfim = kmax + 3*dk

    WRITE(*,*) "ip=",ip
  end do

    Close(2)
    Close(3)
    write(*,*) "k ressonante =", kaux
    write(4,*) "k ressonante =", kaux
    Close(4)
    STOP
  END program wedgearbresv6
  

  SUBROUTINE WALL(x,y,lp,ds,nr)
    REAL(KIND=8) :: x(nr),y(nr),ds,lp
    REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)
    INTEGER :: nr,i

    DO i=1,nr   
      x(i)= 1d0
      y(i)= (DBLE(i)*ds-ds/2.0d0)
    ENDDO
 
    RETURN
  END SUBROUTINE WALL
  
  

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
  
  SUBROUTINE DADOS(R,angcunha,nr,kini,kfim,dk)
    INTEGER nr,aux1,aux2
    REAL(KIND=8) :: R,angcunha,kini,kfim,dk
    CHARACTER(len=6) :: fileid4t
    CHARACTER(len=6) :: fileid5t
  
    aux1 = INT(k)
    aux2 = INT((k-aux1)*1000)
    WRITE(fileid4t,'(i0)') aux1
    WRITE(fileid5t,'(i0)') aux2
    OPEN(unit=10,file='dados'//trim(fileid4t)//'.'//trim(fileid5t)//'.dat',STATUS='UNKNOWN') 
    WRITE(10,*) 'Comprimento', R
    WRITE(10,*) 'angcunha', angcunha
    WRITE(10,*) 'kini', kini
    WRITE(10,*) 'kfim', kfim
    WRITE(10,*) 'dk', dk
    WRITE(10,*) 'nr', nr
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
    END SUBROUTINE MATRIX2
  

  
    SUBROUTINE imagem(x,y,l,angcunha,result) 
    REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)
    REAL(KIND=8) :: x,y,angcunha, ref(2,2),rot(2,2), aux(2),result(2)
    INTEGER      :: l
  
    aux(1) = x
    aux(2) = y
  
    if ( Mod(l,2)==0 ) then
      rot(1,1)= Cos(angcunha*DBLE(l))
      rot(2,2)= rot(1,1)
      rot(1,2)= -Sin(angcunha*DBLE(l))
      rot(2,1)= -rot(1,2)
      result = matmul(rot,aux)
    else
      ref(1,1)= Cos(angcunha*DBLE(l+1))
      ref(2,2)= -ref(1,1)
      ref(1,2)= Sin(angcunha*DBLE(l+1))
      ref(2,1)= ref(1,2)
      result = matmul(ref,aux)  
    end if
    RETURN 
  END SUBROUTINE imagem


  SUBROUTINE submedia(m,nr,media)
    INTEGER :: nr
    REAL(KIND=8) :: media
    COMPLEX(KIND=8) :: m(nr,nr)
    DO i=1,nr
      DO j=1,nr
         media = media + REAL(m(i,j)*DCONJG((m(i,j))))
      ENDDO
   ENDDO
   media = media/(nr**2)
  END SUBROUTINE submedia


  SUBROUTINE subdp(m,nr,media,desviopadrao)
    INTEGER :: nr
    REAL(KIND=8) :: media, desviopadrao
    COMPLEX(KIND=8) :: m(nr,nr)
    
    desviopadrao=0.0d0
    DO i=1,nr
      DO j=1,nr
        desviopadrao =  desviopadrao+( REAL(m(i,j)*DCONJG((m(i,j))))-media)**2
      ENDDO
   ENDDO
   desviopadrao = DSQRT(desviopadrao/((nr**2)-1))
  END SUBROUTINE subdp
