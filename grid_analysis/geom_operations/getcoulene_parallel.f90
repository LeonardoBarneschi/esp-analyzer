      Subroutine getcoulene(atc,atq,protc,protq,nat,nprot,ene)

      IMPLICIT None
      INTEGER, intent(in) :: nat,nprot
      REAL(8), DIMENSION(nat), intent(in) :: atq
      REAL(8), DIMENSION(nat,3), intent(in) :: atc
      REAL(8), DIMENSION(nprot), intent(in) :: protq
      REAL(8), DIMENSION(nprot,3), intent(in) :: protc
      REAL(8), DIMENSION(nat) :: ene_local
      REAL(8), DIMENSION(nat), intent(out) :: ene
      REAL(8), DIMENSION(3) :: rij
      REAL(8) :: r, ang2au
      INTEGER :: i,j,nthrd,omp_get_max_threads

      ang2au=1/0.52917721
      ene=0.d0

!$omp parallel private(i,j,rij,r,ene_local)
      ene_local=0.d0

!$omp do
      do i=1,nat
          do j=1,nprot
              rij=(atc(i,:)-protc(j,:))
              r=dsqrt(rij(1)**2.0+rij(2)**2.0+rij(3)**2.0)
              r=r*ang2au
              ene_local(i)=ene_local(i)+atq(i)*protq(j)/r
          end do
      end do
!$omp end do

!$omp critical
      ene=ene+ene_local
!$omp end critical
!$omp end parallel

      end Subroutine
