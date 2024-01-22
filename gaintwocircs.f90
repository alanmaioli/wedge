!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   Objective: Write in a file the values of the circles parameters, and the respecive C
! mean values related to the segments of the circles.                                  C
!   md is related to the segment of the small circle close to the wedge vertex.        C
!   mi is related to the segments of the two circles.                                  C
!   me is related to the segment of the larger circle, farther from the wedge vertex.  C

program ganhodoiscircsv2
   implicit none
   
REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)
REAL(KIND=8), PARAMETER :: alfa = pi/8.d0
INTEGER, PARAMETER :: nb = 20
REAL(KIND=8), PARAMETER :: R = 0.1d0
REAL(KIND=8), PARAMETER :: razao = 0.05d0
REAL(KIND=8), PARAMETER :: tetac = pi/4.0d0
!----------------------------------------------------------------|
!     Declaração de variáveis                                    |
!----------------------------------------------------------------|
INTEGER i, j, i2, j2, ik, ig, nr, i3, j3, nr1, nr2, aux1, n1, n2
REAL(KIND=8) :: ds, xc1, yc1, xc2, yc2, x0, y0, k, k0, gama0, gama, R2, D
REAL(KIND=8) :: lp, lcx, media, media1, media2, media3, l1, l2


COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: M, mi, me, md

REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xa, ya
CHARACTER(len=6) :: fileid1t

OPEN (3, file='feliz/gama50resmaisalta3.dat', STATUS='UNKNOWN')
write(3,*) "k ","R1 ","xc1  ","yc1  ","R2 ","xc2  ","yc2  ","gama ","nr ","MediaDentro  ","MediaIntermediario ", "MediaExterior"

lcx = 3.0d0
gama0 = 50d0
k0 = 18.5d0
x0 = 0.1d0


do i = 1, 2000
   xc1 = x0 + DBLE(i - 1)*0.00025d0
   yc1 = xc1*DTAN(tetac/2d0)
   D = R 
   WRITE (*, *) "xc1 =", xc1
   if (DABS(DTAN(tetac)*xc1 - yc1)/DSQRT(DTAN(tetac)**2 + 1) .LE. R) then
      write (*, *) "Circulo 1 ultrapassou a parte superior da cunha"
      CYCLE
    
   end if


   xc2 = 1.1d0
   yc2 = xc2*DTAN(tetac/2d0)
   R2 = yc2 - D
   l1 = 2d0*pi*R
   l2 = 2d0*pi*R2
   lp = l1 + l2
  

   if (yc2 > DTAN(tetac)*xc2) then
      CYCLE
   end if

   if (DABS(DTAN(tetac)*xc2 - yc2)/DSQRT(DTAN(tetac)**2 + 1) .LE. R2) then
      CYCLE
   end if

   if (DABS(yc1) .LE. R) then
      CYCLE
   end if

   if (DABS(yc2) .LE. R2) then
      CYCLE
   end if

   if (DSQRT((xc1 - xc2)**2 + (yc1 - yc2)**2) .LE. R2 + R) then
      CYCLE
   end if

   if (xc1 > xc2) then
      CYCLE
   end if

   do ig = 1, 1
      gama = gama0 + dble(ig - 1)*0.02d0

      aux1 = INT(gama)

      do ik = 1, 500
         k = k0 + dble(ik - 1)*0.001d0

         nr = INT(lp*k/(2d0*pi*razao))
         ds = lp/DBLE(nr)
         n1 = int(l1/ds)
         n2 = nr - n1

         ALLOCATE (xa(1:nr), ya(1:nr), m(1:nr, 1:nr), me(1:nr, 1:nr), mi(1:nr, 1:nr), md(1:nr, 1:nr))

         m = 0
         CALL WALL(xa, ya, nr, R, xc1, yc1, R2, xc2, yc2, n1, n2)
         CALL MATRIX2(m, k, ds, xa, ya, nr, tetac, gama)
         CALL identifica(m, nr, mi, md, me, xa, ya, xc1, yc1, xc2, yc2, R, R2)

         CALL submedia2(md, nr, media)
         media1 = media
         CALL submedia2(mi, nr, media)
         media2 = media
         CALL submedia2(me, nr, media)
         media3 = media


         write (3, *) k, R, xc1, yc1, R2, xc2, yc2, gama, nr, media1, media2, media3

         DEALLOCATE (xa, ya, M, me, md, mi)

      end do

   end do
end do

write (*, *) "FIM DO PROGRAMA!!!"

close (3)

end program ganhodoiscircsv2



SUBROUTINE WALL(x, y, nr, R, xc1, yc1, R2, xc2, yc2, n1, n2)
   REAL(KIND=8) :: x(nr), y(nr), R, xc1, yc1, R2, xc2, yc2, inc1, inc2
   REAL(KIND=8), PARAMETER :: pi = DACOS(-1d0)
   INTEGER :: nr, i, n1, n2

   inc1 = 2d0*pi/DBLE(n1)
   inc2 = 2d0*pi/DBLE(n2)
   DO i = 1, nr
      if (i <= n1) then
         x(i) = xc1 + R*DCOS((DBLE(i)*inc1 - inc1/2.0d0))
         y(i) = yc1 + R*DSIN((DBLE(i)*inc1 - inc1/2.0d0))
      else
         x(i) = xc2 + R2*DCOS((DBLE(i - n1)*inc2 - inc2/2.0d0))
         y(i) = yc2 + R2*DSIN((DBLE(i - n1)*inc2 - inc2/2.0d0))
      end if

   END DO

   RETURN
END SUBROUTINE WALL


SUBROUTINE SPACE(xs, ys, lcx, tetac, n)
   REAL(KIND=8) :: xs(n, n), ys(n, n), lcx, lcy, tetac, dsx, dsy
   INTEGER      :: n, i, imax
   lcy = lcx 
   dsx = lcx/DBLE(n)
   dsy = lcy/DBLE(n)
   DO i = 1, n
      imax = min(n, INT(DTAN(tetac)*i)) 
      DO j = 1, imax
         xs(i, j) = DBLE(i)*dsx - dsx/2.0d0
         ys(i, j) = DBLE(j)*dsy - dsy/2.0d0
      END DO
   END DO
END SUBROUTINE SPACE


SUBROUTINE dgreta(ds, k, diag)
   COMPLEX(KIND=8) :: diag, auxc, int
   REAL(KIND=8) :: den, parc, y1, y2, pi, ds, k
   INTEGER :: n
   REAL(KIND=8), PARAMETER :: euler = 0.5772156649015329d0
   pi = DACOS(-1.d0)
   n = 0
   den = 1.d0
   parc = 0.d0
   int = (0.d0, 0.d0)
700 CONTINUE
   den = den*DBLE(2*n)**2
   IF (n == 0) den = 1.d0
   y1 = (2.d0/(2.d0*DBLE(n) + 1.d0))*(ds/2.d0)**(2*n + 1)*k**(2*n)/den
   IF (n /= 0) parc = parc + 1.d0/DBLE(n)
   y2 = (2.d0/pi)*y1*(euler - parc - 1.d0/(2.d0*DBLE(n) + 1.d0) + DLOG(k*ds/4.d0))
   auxc = DCMPLX(y1, y2)
   int = int + (-1.d0)**n*(auxc)
   n = n + 1
   IF (n >= 100) THEN
      WRITE (*, *) 'Muitas iteracoes para a diagonal 1'
      WRITE (*, *) y1, y2, int
      STOP
   END IF
   IF (DABS(y1) > 1.d-20 .OR. DABS(y2) > 1.d-20) GOTO 700
   diag = int
   RETURN
END SUBROUTINE dgreta


SUBROUTINE MATRIX2(m, k, ds, xa, ya, nr, angcunha, gama)

   COMPLEX(KIND=8), PARAMETER :: sigma = (0d0, -0.25d0)
   REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)

   INTEGER :: i, j, nr, ipiv(nr), info, l, Maiusculo
   REAL(KIND=8) :: k, xa(nr), ya(nr), y1, y2, ds, x1, gama
   REAL(KIND=8) :: raux(2), angcunha
   COMPLEX(KIND=8) :: intg, diag, G0, m(nr, nr), work(nr, nr), G1, identidade(nr, nr)

   CALL dgreta(ds, k, diag)

   Maiusculo = Int(pi/angcunha)
   DO i = 1, nr
      DO j = 1, i
         IF (i == j) THEN
            G1 = diag/ds
            G0 = (0.0d0, 0.0d0)
            do l = 1, 2*Maiusculo - 1
               raux = (0.0d0, 0.0d0)
               Call imagem(xa(j), ya(j), l, angcunha, raux)
               x1 = k*DSQRT((raux(1) - xa(i))**2 + (raux(2) - ya(i))**2)
               y1 = bessel_j0(x1)
               y2 = bessel_y0(x1)
               G0 = G0 + ((-1.0d0)**l)*DCMPLX(y1, y2)
            end do
            G0 = G0 + G1


         ELSE

            G0 = (0.0d0, 0.0d0)
            do l = 0, 2*Maiusculo - 1
               raux = (0.0d0, 0.0d0)
               Call imagem(xa(j), ya(j), l, angcunha, raux)

               x1 = k*DSQRT((raux(1) - xa(i))**2 + (raux(2) - ya(i))**2)
               y1 = bessel_j0(x1)
               y2 = bessel_y0(x1)
               G0 = G0 + ((-1.0d0)**l)*DCMPLX(y1, y2)

            end do

         END IF
         intg = sigma*G0*ds
         m(i, j) = intg
         m(j, i) = intg
         IF (CDABS(m(i, j)) == 0.d0) THEN
            WRITE (*, *) i, j, m(i, j), x1
            WRITE (*, *) intg, ds, G0
            stop 'Problemas'
         END IF
      END DO
   END DO


   DO i = 1, nr
      DO j = 1, nr
         if (i /= j) then
            identidade(i, j) = 0
         else
            identidade(i, j) = 1
         end if
      END DO
   END DO

   m = identidade - gama*m


   CALL ZGETRF(nr, nr, m, nr, ipiv, info)
   CALL ZGETRI(nr, m, nr, ipiv, work, nr, info)

   m = gama*m
   RETURN
END SUBROUTINE matrix2



SUBROUTINE imagem(x, y, l, angcunha, result)
   REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)
   REAL(KIND=8) :: x, y, angcunha, ref(2, 2), rot(2, 2), aux(2), result(2)
   INTEGER      :: l

   aux(1) = x
   aux(2) = y

   if (Mod(l, 2) == 0) then
      rot(1, 1) = Cos(angcunha*DBLE(l))
      rot(2, 2) = rot(1, 1)
      rot(1, 2) = -Sin(angcunha*DBLE(l))
      rot(2, 1) = -rot(1, 2)
      result = matmul(rot, aux)
   else
      ref(1, 1) = Cos(angcunha*DBLE(l + 1))
      ref(2, 2) = -ref(1, 1)
      ref(1, 2) = Sin(angcunha*DBLE(l + 1))
      ref(2, 1) = ref(1, 2)
      result = matmul(ref, aux)
   end if
   RETURN
END SUBROUTINE imagem


subroutine identifica(m, nr, mi, md, me, xa, ya, xc1, yc1, xc2, yc2, R, R2)
   REAL(KIND=8), PARAMETER :: pi = DACOS(-1.0d0)
   INTEGER :: nr, i, j
   COMPLEX(KIND=8) :: m(nr, nr), mi(nr, nr), md(nr, nr), me(nr, nr)
   REAL(KIND=8) :: teste1, teste2, aux
   REAL(KIND=8) :: xa(nr), ya(nr), xc1, yc1, xc2, yc2, R, R2

   mi = 0
   md = 0
   me = 0
   do i = 1, nr
      do j = 1, nr

         if (xa(i) <= xc1 + r .AND. xa(j) <= xc1 + r) then

            teste1 = ya(i) + xa(i)*(xc1/yc1) - (xc1**2/yc1) - yc1
            teste2 = ya(j) + xa(j)*(xc1/yc1) - (xc1**2/yc1) - yc1
            if (teste1 < 0d0 .AND. teste2 < 0d0) then
               md(i, j) = m(i, j)
            end if

         else if (xa(i) >= xc2 - R2 .AND. xa(j) >= xc2 - R2) then

            teste1 = ya(i) + xa(i)*(xc2/yc2) - (xc2**2/yc2) - yc2
            teste2 = ya(j) + xa(j)*(xc2/yc2) - (xc2**2/yc2) - yc2
            if (teste1 > 0d0 .AND. teste2 > 0d0) then
               me(i, j) = m(i, j)
            end if
         end if

         aux = (yc2 - yc1)/(xc2 - xc1)
         teste1 = ya(i) + xa(i)/aux - xc1/aux - yc1
         teste2 = ya(i) + xa(i)/aux - xc2/aux - yc2
         if (teste1 >= 0 .AND. teste2 <= 0) then
            teste1 = ya(j) + xa(j)/aux - xc1/aux - yc1
            teste2 = ya(j) + xa(j)/aux - xc2/aux - yc2
            if (teste1 >= 0 .AND. teste2 <= 0) then
               mi(i, j) = m(i, j)
            end if
         end if
      end do
   end do

end subroutine identifica


SUBROUTINE submedia(m, nr, media)
   INTEGER :: nr
   REAL(KIND=8) :: media
   COMPLEX(KIND=8) :: m(nr, nr)
   media = 0.0d0
   DO i = 1, nr
      DO j = 1, nr
         media = media + REAL(m(i, j)*DCONJG((m(i, j))))
      END DO
   END DO
   media = media/(nr**2)
END SUBROUTINE submedia


SUBROUTINE submedia2(maux, nr, media)
   INTEGER :: nr, cont
   REAL(KIND=8) :: media, temp
   COMPLEX(KIND=8) :: maux(nr, nr)
   media = 0.0d0
   cont = 0
   DO i = 1, nr
      DO j = 1, nr
         temp = REAL(maux(i, j)*DCONJG((maux(i, j))))
         if (temp /= 0d0) then
            media = media + temp
            cont = cont + 1
         end if
      END DO
   END DO
   media = media/(DBLE(cont)**2)
END SUBROUTINE submedia2


SUBROUTINE subdp(m, nr, media, desviopadrao)
   INTEGER :: nr
   REAL(KIND=8) :: media, desviopadrao
   COMPLEX(KIND=8) :: m(nr, nr)

   desviopadrao = 0.0d0
   DO i = 1, nr
      DO j = 1, nr
         desviopadrao = desviopadrao + (REAL(m(i, j)*DCONJG((m(i, j)))) - media)**2
      END DO
   END DO
   desviopadrao = DSQRT(desviopadrao/((nr**2) - 1))
END SUBROUTINE subdp
