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
echo "Input the angle of viewpoint:"
read angle

ncl 1> /dev/null <<-EOF
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    numPoint = ${n}
    lonPoint = (/${lon}/)
    latPoint = (/${lat}/)

    lonMin = min(lonPoint)
    lonMax = max(lonPoint)
    latMin = min(latPoint)
    latMax = max(latPoint)

    lonCenter = (lonMin+lonMax)/2
    latCenter = (latMin+latMax)/2

    wks = gsn_open_wks("x11", "spherical_polygon")
    
    mapRes = True
    mapRes@gsnFrame = False
    mapRes@mpGreatCircleLinesOn = True
    mapRes@mpGridAndLimbOn = True
    mapRes@mpGridLineColor = "Background"
    mapRes@mpProjection = "Satellite"
    mapRes@mpCenterLonF = lonCenter
    mapRes@mpCenterLatF = latCenter
    mapRes@mpLimitMode = "Angles"
    mapRes@mpLeftAngleF = ${angle}
    mapRes@mpRightAngleF = ${angle}
    mapRes@mpBottomAngleF = ${angle}
    mapRes@mpTopAngleF = ${angle}

    map = gsn_csm_map(wks, mapRes)

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

    frame(wks)

end
EOF
