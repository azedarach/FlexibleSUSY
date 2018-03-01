program run_loops
  use COLLIER

  Implicit none

  double precision :: p2, m12, m22, mu2
  double complex :: A01,A02,B0,B1,B00
  double complex, allocatable :: Bcoeff(:,:),Bcoeffuv(:,:)
  double precision, allocatable :: Berr(:)
  integer :: i, rank
  character(len=255) :: read_line, file_name
  logical :: file_exists=.False.

  ! for B coefficients
  rank = 3
  allocate(Bcoeff(0:rank/2,0:rank))
  allocate(Bcoeffuv(0:rank/2,0:rank))
  allocate(Berr(0:rank))

  ! open input file
  call get_command_argument(1, file_name)
  inquire(file=trim(file_name), exist=file_exists)
  if (.not.file_exists) then
     write(*,*) "Error: file ",trim(file_name)," does not exist!"
     stop 1
  end if
  open(1, file=trim(file_name), err=1)

  call Init_cll(6,6,"")

  write (*,*) "# A0(m0^2) A0(m1^2) B0(p^2,m0^2,m1^2) B1(p^2,m0^2,m1^2) B00(p^2,m0^2,m1^2) B11(p^2,m0^2,m1^2)"

  do
     read(1,"(a255)", end=1, err=1) read_line
     if (read_line(1:1).eq."#") cycle
     if (read_line.eq." ") cycle

     read(1,*) p2, m12, m22, mu2

     call SetMuUV2_cll(mu2)
     call SetMuIR2_cll(mu2)

     call A0_cll(A01,dcmplx(m12))
     call A0_cll(A02,dcmplx(m22))
     call B_cll(Bcoeff,Bcoeffuv,dcmplx(p2),dcmplx(m12),dcmplx(m22),rank,Berr)  

     B0 = Bcoeff(0,0)
     B1 = Bcoeff(0,1)
     B00 = Bcoeff(1,0)

     write(*,*) real(A01), real(A02), real(B0), real(B1), real(B00)
  end do

1 close(1)

end program run_loops
