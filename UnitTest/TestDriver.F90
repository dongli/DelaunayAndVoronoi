program TestDriver

    use SampleManager
    use SphereService
    use DelaunayAndVoronoi

    implicit none

    integer, parameter :: n = 9
    real(8) lon(n), lat(n)
    real(8) dlon
    integer i

    dlon = PI2/(n-1)
    lon = 0.0d0
    lat = PI05
    do i = 2, n
        lon(i) = (i-2)*dlon
        lat(i) = PI/3.0d0
    end do

    call SampleManager_Init(n)
    call SampleManager_Set(lon, lat)
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation
    call DelaunayAndVoronoi_Output("test.nc")

end program TestDriver
