#!/bin/bash

echo "Input the coordinates of the vertices:"
lon1=go
n=-1
while [ $lon1 ]; do
    read lon1 lat1
    lon=$lon,$lon1
    lat=$lat,$lat1
    n=$[$n+1]
done
lon=${lon%%,}
lon=${lon##,}
lat=${lat%%,}
lat=${lat##,}

mapProj="Satellite"
echo "Input map projection:"
read -p "[Default: $mapProj]> " ans
if [ -n "$ans" ]; then
    mapProj=$ans
fi

lonCenter=180
latCenter=90
angle="1.0"
if [ "$mapProj" == "Satellite" ]; then
    echo "Input the center of viewpoint:"
    read -p "> " lonCenter latCenter
    echo "Input the angle of viewpoint:"
    read -p "[Default:$angle]> " ans
    if [ -n "$ans" ]; then
        angle=$ans
    fi
fi

ncl <<-EOF
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    numPoint = $n
    lonPoint = (/$lon/)
    latPoint = (/$lat/)

    lonMin = min(lonPoint)
    lonMax = max(lonPoint)
    latMin = min(latPoint)
    latMax = max(latPoint)

    lonCenter = (lonMin+lonMax)/2
    latCenter = (latMin+latMax)/2

    wks = gsn_open_wks("pdf", "spherical_polygon")
  
    mapProj = "$mapProj"
  
    mapRes = True
    mapRes@gsnDraw = False
    mapRes@gsnFrame = False
    mapRes@mpGreatCircleLinesOn = True
    mapRes@mpGridAndLimbOn = True
    mapRes@mpGridLineColor = "Background"

    if (mapProj .eq. "Satellite") then
        mapRes@mpProjection = "Satellite"
        mapRes@mpCenterLonF = $lonCenter
        mapRes@mpCenterLatF = $latCenter
        mapRes@mpLimitMode = "Angles"
        mapRes@mpLeftAngleF = $angle
        mapRes@mpRightAngleF = $angle
        mapRes@mpBottomAngleF = $angle
        mapRes@mpTopAngleF = $angle
    end if

    map = gsn_csm_map(wks, mapRes)

    textRes = True
    textRes@txFontHeightF = 0.01
    textRes@txFont = "helvetica-bold"

    do i = 0, numPoint-1
        text = gsn_add_text(wks, map, sprinti("%d", i+1), lonPoint(i)-0.5, latPoint(i)-0.5, textRes)
    end do

    draw(map)

    edgeRes = True
    edgeRes@gsLineThicknessF = 0.5
    edgeRes@gsLineColor = "blue"

    lon = new(2, float)
    lat = new(2, float)

    do i = 0, numPoint-1
        j = i-1
        if (i .eq. 0) then
            j = numPoint-1
        end if
        lon(0) = lonPoint(j)
        lon(1) = lonPoint(i)
        lat(0) = latPoint(j)
        lat(1) = latPoint(i)
        gsn_polyline(wks, map, lon, lat, edgeRes)
    end do

    vertexRes = True
    vertexRes@gsMarkerIndex = 1
    vertexRes@gsMarkerSizeF = 0.01
    vertexRes@gsMarkerColor = "red"

    gsn_polymarker(wks, map, lonPoint, latPoint, vertexRes)

    if (False) then
        fillRes = True
        fillRes@gsFillColor = "red"

        x = new(numPoint, float)
        y = new(numPoint, float)
        datatondc(map, lonPoint, latPoint, x, y)
        gsn_polygon_ndc(wks, x, y, fillRes)
    end if

    frame(wks)

end
EOF
