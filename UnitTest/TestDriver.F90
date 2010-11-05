program TestDriver

    use SampleManager
    use SphereService
    use DelaunayAndVoronoi

    implicit none

    integer, parameter :: n = 9
    real(8) lon(n), lat(n)
    real(8) dlon
    integer i

    if (.false.) then
        ! Case 1
        dlon = PI2/(n-1)
        lon = 0.0d0
        lat = PI05
        do i = 2, n
            lon(i) = (i-2)*dlon
            lat(i) = PI/3.0d0
        end do
    else
        ! Case 2
        lon = [45.0,30.0,60.0,45.0,60.0,70.0,31.0,45.0,80.0]/Rad2Deg
        lat = [60.0,30.0,30.0,45.0,10.0,55.0,40.0,50.0,35.0]/Rad2Deg
    end if

    call SampleManager_Init(n)
    call SampleManager_Set(lon, lat)
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation
    call ExtractVoronoiDiagram
    !call DelaunayAndVoronoi_Report
    call DelaunayAndVoronoi_Output("test.nc")

end program TestDriver
