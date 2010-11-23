program TestDriver

    use SampleManager
    use SphereService
    use DelaunayAndVoronoi

    implicit none

    integer n, nx, ny
    real(RealKind), allocatable :: lon(:), lat(:)
    real(RealKind) dlon, dlat, dlon05, dlat05
    integer i, j, k

    if (.false.) then
        ! Case 1
        n = 9
        allocate(lon(n))
        allocate(lat(n))
        dlon = PI2/(n-1)
        lon = 0.0d0
        lat = PI05
        do i = 2, n
            lon(i) = (i-2)*dlon
            lat(i) = PI/3.0d0
        end do
    else if (.false.) then
        ! Case 2
        n = 9
        allocate(lon(n))
        allocate(lat(n))
        lon = [45.0,30.0,60.0,45.0,60.0,70.0,31.0,45.0,80.0]/Rad2Deg
        lat = [60.0,30.0,30.0,45.0,10.0,55.0,40.0,50.0,35.0]/Rad2Deg
    else
        nx = 100
        ny = 100
        n = nx*ny
        allocate(lon(n))
        allocate(lat(n))
        dlon = PI2/nx
        dlat = PI/ny
        dlon05 = dlon*0.5
        dlat05 = dlat*0.5
        k = 1
        do j = 1, ny
            do i = 1, nx
                lon(k) = dlon05+(i-1)*dlon
                lat(k) = PI05-dlat05-(j-1)*dlat
                k = k+1
            end do
        end do
    end if

    call SampleManager_Init(n)
    call SampleManager_Set(lon, lat)
    !doRandom = .false.
    !idx1 = 19; idx2 = 20; idx3 = 75
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation
    call ExtractVoronoiDiagram
#if (defined KeepVirtualDVT)
    call DelaunayAndVoronoi_OutputDelaunay("only_DT_output.nc")
#else
    call DelaunayAndVoronoi_Output("output.nc")
#endif

end program TestDriver
