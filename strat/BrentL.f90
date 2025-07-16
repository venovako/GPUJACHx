! ifort.exe /standard-semantics /fast BrentL.f90 /link /RELEASE
! ifort -standard-semantics -fast BrentL.f90 -o BrentL.exe -s

program BrentL

  implicit none

  character(len=6) :: arg1
  integer :: n, n1, n2, n_2, k, s, i, j
  integer, allocatable :: L(:), R(:), out_L(:), out_R(:), in_L(:), in_R(:)

  integer, intrinsic :: mod, command_argument_count
  intrinsic :: get_command_argument

  if (command_argument_count() .ne. 1) stop 'BrentL.exe n'
  call get_command_argument(1,arg1)
  read (arg1,*) n

  if (n .lt. 2) stop 'n < 2'
  n1 = n - 1
  n2 = n - 2
  if (mod(n,2) .ne. 0) stop 'n odd'
  n_2 = n / 2

  allocate(L(n_2))
  allocate(R(n_2))
  allocate(out_L(n_2))
  allocate(out_R(n_2))
  allocate(in_L(n_2))
  allocate(in_R(n_2))

  write (*,1) '#ifdef USE_STRAT_ARRAY_DECLARATOR'
  write (*,1,advance='no') 'unsigned '
  if (n .le. 256) then
     write (*,1,advance='no') 'char  '
  else
     write (*,1,advance='no') 'short '
  end if
  write (*,2) n, n1, n_2
  write (*,1) '#endif /* USE_STRAT_ARRAY_DECLARATOR */'
  write (*,1) '{'

  ! init
  do k = 1, n_2
     L(k) = 2 * k - 1
     R(k) = 2 * k
  end do

  ! unused
  out_L(1) = 0
  in_L(1) = 0
  out_R(n_2) = 0
  in_R(n_2) = 0

  do s = 0, n2
     write (*,1,advance='no') '  {'
     do k = 1, n_2
        if (L(k) .lt. R(k)) then
           i = L(k) - 1
           j = R(k) - 1
        else
           i = R(k) - 1
           j = L(k) - 1
        end if
        if (n .le. 10) then
           write (*,3,advance='no') i, j
        else if (n .le. 100) then
           write (*,4,advance='no') i, j
        else if (n .le. 1000) then
           write (*,5,advance='no') i, j
        else if (n .le. 10000) then
           write (*,6,advance='no') i, j
        else
           write (*,7,advance='no') i, j
        end if
        if (k .lt. n_2) write (*,1,advance='no') ','
        if (k .eq. 1) then
           out_R(k) = R(k)
        else if (k .lt. n_2) then
           out_R(k) = L(k)
        end if
        if (k .gt. 1) out_L(k) = R(k)
     end do
     do k = 1, n_2 - 1
        in_L(k + 1) = out_R(k)
     end do
     do k = 1, n_2 - 1
        in_R(k) = out_L(k + 1)
     end do
     do k = 1, n_2
        if (k .lt. n_2) then
           R(k) = in_R(k)
        else
           R(k) = L(k)
        end if
        if (k .gt. 1) L(k) = in_L(k)
     end do
     if (s .lt. n2) then
        write (*,1) '},'
     else
        write (*,1) '}'
     end if
  end do

  write (*,1) '};'

  deallocate(in_R)
  deallocate(in_L)
  deallocate(out_R)
  deallocate(out_L)
  deallocate(R)
  deallocate(L)

1 format(a)
2 format('BrentL',i5.5,'[',i5,'][',i5,'][2] =')
3 format('{',i1,',',i1,'}')
4 format('{',i2,',',i2,'}')
5 format('{',i3,',',i3,'}')
6 format('{',i4,',',i4,'}')
7 format('{',i5,',',i5,'}')

end program BrentL
