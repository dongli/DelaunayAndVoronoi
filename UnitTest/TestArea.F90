program TestArea

    use SphereService

    implicit none

    integer, parameter :: n = 4
    real(RealKind) lon(n), lat(n)
    real(RealKind) x(n), y(n), z(n)
    real(RealKind) area

    integer i

    lon = [  0.0,  0.0,  0.0, 90.0]/Rad2Deg
    lat = [ 90.0,  0.0,-90.0,  0.0]/Rad2Deg

    do i = 1, n
        call CartesianTransformOnUnitSphere(lon(i), lat(i), &
            x(i), y(i), z(i))
    end do

    write(*, "('The area of a spherical polygon:')")
    call calcArea(n, x, y, z, area)
    write(*, *) area/PI/Re2, "PI"
    call calcArea(n, lon, lat, area)
    write(*, *) area/PI/Re2, "PI"

end program
