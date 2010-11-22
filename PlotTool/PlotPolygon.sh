#!/bin/bash

which ncl > /dev/null
if [ $? != 0 ]; then
    echo "*** You need to install NCL to run this script!"
    exit -1
fi

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

mapProj="CE"
echo "Input map projection (CE/ST):"
read -p "[Default: $mapProj] > " ans
if [ -n "$ans" ]; then
    mapProj=$ans
fi

lonCnt="180"
latCnt="0"
angle="1"

if [ $mapProj == "CE" ]; then
    echo "** Cylindrical equidistant projection is used."
fi

if [ $mapProj == "ST" ]; then
    echo "** Satellite projection is used."
    echo "Input the center of viewport:"
    read -p "[Default: lon=$lonCnt, lat=$latCnt] > " ans1 ans2
    if [ -n "$ans1" -a -n "$ans2" ]; then
        lonCnt=$ans1
        latCnt=$ans2
    fi
    echo "Input the angle of viewport:"
    read -p "[Default: $angle] > " ans
    if [ -n "$ans" ]; then
        angle=$ans
    fi
fi

figtype="x11"
echo "Input figure type: (x11/pdf)"
read -p "[Default: $figtype] > " ans
if [ -n "$ans" ]; then
    figtype=$ans
fi

ncl <<-EOF
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    numPnt = $n
    lonPnt = (/$lon/)
    latPnt = (/$lat/)

    lonMin = min(lonPnt)
    lonMax = max(lonPnt)
    latMin = min(latPnt)
    latMax = max(latPnt)

    lonCnt = (lonMin+lonMax)/2
    latCnt = (latMin+latMax)/2

    wks = gsn_open_wks("$figtype", "spherical_polygon")
  
    mapProj = "$mapProj"
  
    mapRes = True
    mapRes@gsnDraw = False
    mapRes@gsnFrame = False
    mapRes@mpGreatCircleLinesOn = True
    mapRes@mpGridAndLimbOn = True
    mapRes@mpGridLineColor = "Background"

    if (mapProj .eq. "ST") then
        mapRes@mpProjection = "Satellite"
        mapRes@mpCenterLonF = $lonCnt
        mapRes@mpCenterLatF = $latCnt
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

    do i = 0, numPnt-1
        text = gsn_add_text(wks, map, sprinti("%d", i+1), lonPnt(i)-0.5, latPnt(i)-0.5, textRes)
    end do

    draw(map)

    edgeRes = True
    edgeRes@gsLineThicknessF = 0.5
    edgeRes@gsLineColor = "blue"

    lon = new(2, float)
    lat = new(2, float)

    do i = 0, numPnt-1
        j = i-1
        if (i .eq. 0) then
            j = numPnt-1
        end if
        lon(0) = lonPnt(j)
        lon(1) = lonPnt(i)
        lat(0) = latPnt(j)
        lat(1) = latPnt(i)
        gsn_polyline(wks, map, lon, lat, edgeRes)
    end do

    vertexRes = True
    vertexRes@gsMarkerIndex = 1
    vertexRes@gsMarkerSizeF = 0.01
    vertexRes@gsMarkerColor = "red"

    gsn_polymarker(wks, map, lonPnt, latPnt, vertexRes)

    if (False) then
        fillRes = True
        fillRes@gsFillColor = "red"

        x = new(numPnt, float)
        y = new(numPnt, float)
        datatondc(map, lonPnt, latPnt, x, y)
        gsn_polygon_ndc(wks, x, y, fillRes)
    end if

    frame(wks)

end
EOF
