program CompareWithSTRIPACK

    use SphereService
    use SampleManager
    use DelaunayAndVoronoi

    implicit none

    integer, parameter :: nx = 10, ny = 10, n = 9!nx*ny
    real(RealKind) lon(n), lat(n)
    real(RealKind) dlon, dlat, dlon05, dlat05
    integer i, j, k

    ! Internal variables for STRIPACK
    integer, parameter :: n6 = 6*n
    integer, parameter :: n2 = 2*n
    integer LIST(n6), LPTR(n6), LEND(n), LNEW, NEAR(n), NEXT(n), DIST(n), IER
    real X(n), Y(n), Z(n)
    real ELAT, VLAT, ELON, VLON

    real startTime1, finishTime1 ! CPU timer for DelaunayAndVoronoi
    real startTime2, finishTime2 ! CPU timer for STRIPACK

    if (.false.) then
        dlon = PI2/nx
        dlat = PI/ny
        dlon05 = dlon*0.5
        dlat05 = dlat*0.5
        k = 1
        do j = 1, ny
            do i = 1, nx
                lon(k) = dlon05+i*dlon
                lat(k) = PI05-dlat05-(j-1)*dlat
                k = k+1
            end do
        end do
    else
        ! Case 2
        lon = [45.0,30.0,60.0,45.0,60.0,70.0,31.0,45.0,80.0]/Rad2Deg
        lat = [60.0,30.0,30.0,45.0,10.0,55.0,40.0,50.0,35.0]/Rad2Deg
    end if

    ! Call DelaunayAndVoronoi ........................................... STEP1

    call cpu_time(startTime1)

    call SampleManager_Init(n)
    call SampleManager_Set(lon, lat)
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation

    call cpu_time(finishTime1)

    write(*, "('The used CPU time by DelaunayAndVoronoi: ', F10.6)") &
        finishTime1-startTime1

    call DelaunayAndVoronoi_OutputDelaunay("DelaunayAndVoronoi_output.nc")

    ! Call STRIPACK ..................................................... STEP2

    call cpu_time(startTime2)

    call TRANS(n, lat, lon, X, Y, Z)
    call TRMESH(n, X, Y, Z, LIST, LPTR, LEND, LNEW, NEAR, NEXT, DIST, IER)
    if (IER == -2) then
        write(*, "('STRIPACK: The first three nodes are collinear.')")
        stop
    else if (IER > 0) then
        write(*, "('STRIPACK: Duplicate nodes encountered.')")
        stop
    end if

    call cpu_time(finishTime2)

    write(*, "('The used CPU time by STRIPACK:           ', F10.6)") &
        finishTime2-startTime2

    open(3, file="STRIPACK_output.ps")
    call TRPLOT(3, 7.5, 45.0, 60.0, 90.0, n, X, Y, Z, LIST, LPTR, LEND, &
        "(STRIPACK OUTPUT)", .false., IER)
    close(11)

    stop

end program CompareWithSTRIPACK
