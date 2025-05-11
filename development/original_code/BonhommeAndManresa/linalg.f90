module linalg
  implicit none
contains

  subroutine ludcmp(a,indx,d)
    use utils, only : assert_eq, imaxloc, nrerror, outerprod, swap
    implicit none

    real(8), dimension(:,:), intent(inout) :: a
    integer, dimension(:), intent(out) :: indx
    real(8), intent(out) :: d
    real(8), dimension(size(a,1)) :: vv
    real(8), parameter :: tiny = 1.e-20
    integer :: j,n,imax

    n = assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d = 1.d0
    vv = maxval(abs(a),dim=2)
    if (any(vv == 0.d0)) call nrerror('singular matrix in ludcmp')
    vv = 1.d0 / vv
    do j = 1, n
       imax = (j-1) + imaxloc(vv(j:n)*abs(a(j:n,j)))
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d = -d
          vv(imax) = vv(j)
       end if
       indx(j) = imax
       if (a(j,j) == 0.d0) a(j,j) = tiny
       a(j+1:n,j) = a(j+1:n,j) / a(j,j)
       a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - outerprod(a(j+1:n,j),a(j,j+1:n))
    end do
  end subroutine ludcmp


  subroutine lubksb(a,indx,b)
    use utils, only : assert_eq
    implicit none
    real(8), dimension(:,:), intent(in) :: a
    integer, dimension(:), intent(in) :: indx
    real(8), dimension(:), intent(inout) :: b
    integer :: i, n, ii, ll
    real(8) :: summ

    n = assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii = 0
    do i = 1, n
       ll = indx(i)
       summ = b(ll)
       b(ll) = b(i)
       if (ii /= 0) then
          summ = summ - dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.d0) then
          ii = i
       end if
       b(i) = summ
    end do
    do i = n,1,-1
       b(i) = (b(i) - dot_product(a(i,i+1:n),b(i+1:n))) / a(i,i)
    end do
  end subroutine lubksb


  subroutine lu_inverse(mat,imat)
    use utils, only : unit_matrix, assert_eq
    implicit none
    real(8), dimension(:,:), intent(in) :: mat
    real(8), dimension(:,:), intent(out) :: imat
    real(8), dimension(size(mat,1),size(mat,2)) :: cpmat
    integer, dimension(size(mat,1)) :: indx
    real(8) :: d
    integer :: j,n

    n = assert_eq(size(mat,1),size(mat,2),size(imat,1),size(imat,2),'lu_inverse')
    call unit_matrix(imat,n)
    cpmat = mat
    call ludcmp(cpmat,indx,d)
    do j = 1, n
       call lubksb(cpmat,indx,imat(:,j))
    end do

  end subroutine lu_inverse

end module linalg
