program TestRealData

    use NFWrap
    use SampleManager
    use DelaunayAndVoronoi

    implicit none

    type(FileCard) fcard

    integer numPoint
    real(RealKind), allocatable :: lon(:), lat(:)

    call NFWrap_OpenForRead("C02562.global.nc", fcard)
    call NFWrap_GetDimSize(fcard, "grid_size", numPoint)
    allocate(lon(numPoint), lat(numPoint))
    call NFWrap_Input1DVar(fcard, "grid_center_lon", lon)
    call NFWrap_Input1DVar(fcard, "grid_center_lat", lat)
    call NFWrap_Close(fcard)

    call SampleManager_Init(numPoint)
    call SampleManager_Set(lon, lat)
    call DelaunayAndVoronoi_LinkSample
    call ConstructDelaunayTriangulation
    call ExtractVoronoiDiagram
    call DelaunayAndVoronoi_Output("output.nc")

end program TestRealData
