load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    f = addfile(fileName, "r")

    dims = dimsizes(f->triangle)
    numTriangle = dims(0)

    wks = gsn_open_wks("pdf", figureName)

    ; Plot map
    mapRes = True
    mapRes@gsnFrame = False
    mapRes@mpGreatCircleLinesOn = True
    ;mapRes@mpProjection = "Setereographic"
    mapRes@mpProjection = "Satellite"
    mapRes@mpCenterLonF = 70.0
    mapRes@mpCenterLatF = 60.0

    map = gsn_csm_map(wks, mapRes)

    ; Plot edges
    edgeRes = True
    edgeRes@gsLineColor = "blue"

    lon = new(3, "float")
    lat = new(3, "float")
    do i = 0, numTriangle-1
        do j = 0, 2
            k = f->triangle(i,j)-1
            lon(j) = f->point_lon(k)
            lat(j) = f->point_lat(k)
        end do
        gsn_polyline(wks, map, lon, lat, edgeRes)
    end do

    ; Plot circumcirlces
    circRes = True
    circRes@gsLineThicknessF = 2.
    circRes@gsLineColor = "green"

    arclon = new(50, float)
    arclat = new(50, float)
    do i = 0, numTriangle-1
        nggcog(f->clat(i), f->clon(i), f->radius(i), arclat, arclon)
        gsn_polyline(wks, map, arclon, arclat, circRes)
    end do

    ; Plot vertices
    vertexRes = True
    vertexRes@gsMarkerIndex = 1
    vertexRes@gsMarkerSizeF = 0.02
    vertexRes@gsMarkerColor = "red"

    gsn_polymarker(wks, map, f->point_lon, f->point_lat, vertexRes)

    frame(wks)

end

