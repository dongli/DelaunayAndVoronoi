program TestRealData

    use NFWrap
    use SampleManager
    use DelaunayAndVoronoi

    implicit none

    character(255) filePath
    type(FileCard) fcard

    integer numPoint, numTimeStep
    real(8), allocatable :: lon(:), lat(:)

    call get_command_argument(1, filePath)

    call NFWrap_OpenForRead(filePath, fcard)
    call NFWrap_GetDimSize(fcard, "q_num", numPoint)
    call NFWrap_GetDimSize(fcard, "time", numTimeStep)
    allocate(lon(numPoint), lat(numPoint))
    call NFWrap_Input1DVar(fcard, "q_lon", lon, numTimeStep)
    call NFWrap_Input1DVar(fcard, "q_lat", lat, numTimeStep)

    call SampleManager_Init(numPoint)
    call SampleManager_Set(lon, lat)
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation
    call ExtractVoronoiDiagram
    call DelaunayAndVoronoi_Output("output.nc")

end program TestRealData
