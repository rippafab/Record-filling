program quantile_mapping
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: n, ny, n2, i
    real(dp), allocatable :: time_x(:), x_val(:)
    real(dp), allocatable :: time_y(:), y_val(:)
    real(dp), allocatable :: t1(:), x1(:), m1(:), y1(:)
    
    ! --- 1. Lettura parametri e allocazione ---
    ! Assumiamo di conoscere N o di contarli. Qui contiamo le righe.
    n = count_lines('c1.dat')
    ny = count_lines('c2.dat')
    n2 = count_lines('c3.dat')

    allocate(time_x(n), x_val(n), time_y(ny), y_val(ny))
    allocate(t1(n2), x1(n2), m1(n2), y1(n2))

    ! --- 2. Caricamento dati ---
    open(10, file='c1.dat', status='old')
    do i = 1, n
        read(10, *) time_x(i), x_val(i)
    end do
    close(10)

    open(11, file='c2.dat', status='old')
    do i = 1, ny
        read(11, *) time_y(i), y_val(i)
    end do
    close(11)

    open(12, file='c3.dat', status='old')
    do i = 1, n2
        read(12, *) t1(i), x1(i), m1(i)
    end do
    close(12)

    ! --- 3. Ordinamento delle distribuzioni base ---
    ! Ordiniamo x_val e y_val per costruire le CDF empiriche
    call sort(x_val)
    call sort(y_val)

    ! --- 4. Trasformazione (Interpolazione sui quantili) ---
    ! Mappiamo ogni x1(i) su y1(i) usando il rango in x_val
    do i = 1, n2
        y1(i) = map_value(x1(i), x_val, n, y_val, ny)
    end do

    ! --- 5. Scrittura file finale ---
    open(13, file='c4.dat', status='replace')
    do i = 1, n2
        write(13, '(3F15.6)') t1(i), y1(i), m1(i)
    end do
    close(13)

    print *, "Trasformazione completata. File c4.dat generato."

contains

    ! Funzione per contare le righe di un file
    integer function count_lines(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, ios
        character(len=1024) :: line
        count_lines = 0
        open(newunit=unit, file=filename, status='old')
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) count_lines = count_lines + 1
        end do
        close(unit)
    end function count_lines

    ! Algoritmo di ordinamento semplice (Shell Sort)
    subroutine sort(a)
        real(dp), intent(inout) :: a(:)
        integer :: i, j, h, nn
        real(dp) :: v
        nn = size(a)
        h = 1
        do while (h <= nn/3)
            h = 3*h + 1
        end do
        do while (h >= 1)
            do i = h + 1, nn
                v = a(i)
                j = i
                do while (j > h .and. a(j-h) > v)
                    a(j) = a(j-h)
                    j = j - h
                end do
                a(j) = v
            end do
            h = h / 3
        end do
    end subroutine sort

    ! Funzione di mappatura tramite interpolazione lineare sui ranghi
    real(dp) function map_value(val, sorted_x, nx, sorted_y, ny)
        real(dp), intent(in) :: val
        real(dp), intent(in) :: sorted_x(nx), sorted_y(ny)
        integer, intent(in) :: nx, ny
        integer :: idx, idy
        real(dp) :: frac, p, yfrac, yy

        if (nx <= 1) then
            map_value = sorted_y(max(1, min(ny, 1)))
            return
        end if

        if (ny <= 1) then
            map_value = sorted_y(1)
            return
        end if

        ! Se il valore è fuori range, facciamo clipping
        if (val <= sorted_x(1)) then
            p = 0.0_dp
        else if (val >= sorted_x(nx)) then
            p = 1.0_dp
        else
            idx = binary_search(val, sorted_x, nx)
            if (sorted_x(idx+1) > sorted_x(idx)) then
                frac = (val - sorted_x(idx)) / (sorted_x(idx+1) - sorted_x(idx))
            else
                frac = 0.0_dp
            end if
            p = (real(idx-1, dp) + frac) / real(nx-1, dp)
        end if

        yy = 1.0_dp + p * real(ny-1, dp)
        idy = int(floor(yy))
        if (idy < 1) idy = 1
        if (idy >= ny) then
            map_value = sorted_y(ny)
            return
        end if
        yfrac = yy - real(idy, dp)
        map_value = sorted_y(idy) + yfrac * (sorted_y(idy+1) - sorted_y(idy))
    end function map_value

    integer function binary_search(val, arr, n)
        real(dp), intent(in) :: val, arr(n)
        integer, intent(in) :: n
        integer :: low, high, mid
        low = 1
        high = n
        do while (high - low > 1)
            mid = (low + high) / 2
            if (arr(mid) <= val) then
                low = mid
            else
                high = mid
            end if
        end do
        binary_search = low
    end function binary_search

end program quantile_mapping
