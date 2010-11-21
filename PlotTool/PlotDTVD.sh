#!/bin/bash

which ncl > /dev/null
if [ $? != 0 ]; then
    echo "*** You need to install NCL to run this script!"
    exit -1
fi

ls *.nc > /dev/null 2>&1
if [ $? != 0 ]; then
    echo "*** There is no NetCDF data file!"
    exit -2
fi

echo "Input the NetCDF data file or wildcard:"
echo "Possible files:"
files=$(ls -t *.nc)
echo $files
fileWildCard=$(echo $files | cut -d ' ' -f 1)
read -p "[Default: $fileWildCard] > " ans
if [ -n "$ans" ]; then
    fileWildCard=$ans
fi

DT="True"
echo "Plot Delaunay triangulation (DT)? (y/n)"
read -p "[Default: y] > " ans
if [ -n "$ans" ]; then
    if [ $ans == "n" ]; then
        DT=False
    fi
fi

VD="True"
echo "Plot Voronoi diagram (VD)? (y/n)"
read -p "[Default: y] > " ans
if [ -n "$ans" ]; then
    if [ "$ans" == "n" ]; then
        VD=False
    fi
fi

DVT="True"
echo "Plot vertices? (y/n)"
read -p "[Default: y] > " ans
if [ -n "$ans" ]; then
    if [ "$ans" == "n" ]; then
        DVT=False
    fi
fi

mapProj="CE"
echo "Input map projection (NH/SH/CE/ST):"
read -p "[Default: $mapProj] > " ans
if [ -n "$ans" ]; then
    mapProj=$ans
fi

if [ $mapProj == "CE" ]; then
    echo "** Cylindrical equidistant projection is used."
fi

minLat="80"
if [ $mapProj == "NH" ]; then
    echo "Northen stereographic projection is used."
    echo "Input the southest latitude to view:"
    read -p "[Default: $minLat] > " ans
    if [ -n "$ans" ]; then
        minLat=$ans
    fi
fi

maxLat="-80"
if [ $mapProj == "SH" ]; then
    echo "** Southen stereographic projection is used."
    echo "Input the northest latitude to view:"
    read -p "[Default: $maxLat] > " ans
    if [ -n "$ans" ]; then
        maxLat=$ans
    fi
fi

lonCnt="180"
latCnt="0"
if [ $mapProj == "ST" ]; then
    echo "** Satellite projection is used."
    echo "Input the center of viewport:"
    read -p "[Default: lon=$lonCnt, lat=$latCnt] > " ans1 ans2
    if [ -n "$ans1" -a -n "$ans2" ]; then
        lonCnt=$ans1
        latCnt=$ans2
    fi
fi

figtype="pdf"
echo "Input figure type: (pdf/x11)"
read -p "[Default: $figtype] > " ans
if [ -n "$ans" ]; then
    figtype=$ans
fi

for file in $(ls ${fileWildCard})
do
    figure=$(basename $file .nc)
    echo "Plotting $file ... "
    ncl <<-EOF
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

    begin

        Rad2Deg = 45./atan(1.)

        f = addfile("$file", "r")

        numPnt = dimsizes(f->lonPnt)

        hasDelaunay = False
        if (isfilevar(f, "triIdx")) then
            hasDelaunay = True
        end if

        if (.not. hasDelaunay .and. $DT) then
            system("echo '*** You choose to plot DT, but there is no data.'")
        end if

        lonPnt = f->lonPnt*Rad2Deg
        latPnt = f->latPnt*Rad2Deg

        if (hasDelaunay) then
            numTriangle = dimsizes(f->triIdx(:,0))
        end if

        hasVoronoi = False
        if (isfilevar(f, "numVtx")) then
            hasVoronoi = True
        end if

        if (hasVoronoi) then
            lonVtx = f->lonVtx*Rad2Deg
            latVtx = f->latVtx*Rad2Deg
        end if

        if (.not. hasVoronoi .and. $VD) then
            system("echo '*** You choose to plot VD, but there is no data.'")
        end if

        wks = gsn_open_wks("$figtype", "$figure")

        ; Plot map
        mapRes = True
        mapRes@gsnDraw = False
        mapRes@gsnFrame = False
        mapRes@mpGreatCircleLinesOn = True
        mapRes@mpGridAndLimbOn = True
        mapRes@mpGridLineColor = "Background"

        mapProj = "$mapProj"

        if (mapProj .eq. "SH") then
            mapRes@mpProjection = "Stereographic"
            mapRes@gsnPolar = "SH"
            mapRes@mpMaxLatF = $maxLat
        end if
        if (mapProj .eq. "NH") then
            mapRes@mpProjection = "Stereographic"
            mapRes@gsnPolar = "NH"
            mapRes@mpMinLatF = $minLat
        end if
        if (mapProj .eq. "ST") then
            mapRes@mpProjection = "Satellite"
            mapRes@mpCenterLonF = $lonCnt
            mapRes@mpCenterLatF = $latCnt
        end if

        map = gsn_csm_map(wks, mapRes)

        if (False) then
            ; Plot vertex indices
            textRes = True
            textRes@txFontHeightF = 0.01
            textRes@txFont = "helvetica-bold"

            do i = 0, numPnt-1
                text = gsn_add_text(wks, map, sprinti("%d", i+1), \
                    lonPnt(i)-0.5, latPnt(i)-0.5, textRes)
            end do
        end if

        draw(map)

        if ($DT .and. hasDelaunay) then
            ; Plot Delaunay triangle edges
            system("echo '  Plotting triangle edges ...'")
            edgeRes = True
            edgeRes@gsLineThicknessF = 0.5
            edgeRes@gsLineColor = "blue"

            lon = new(4, "float")
            lat = new(4, "float")
            do i = 0, numTriangle-1
                n = 0
                do j = 0, 2
                    k = f->triIdx(i,j)-1
                    lon(n) = lonPnt(k)
                    lat(n) = latPnt(k)
                    n = n+1
                end do
                if (n .gt. 1) then
                    lon(n) = lon(0)
                    lat(n) = lat(0)
                    gsn_polyline(wks, map, lon(0:n), lat(0:n), edgeRes)
                end if
            end do
            delete(lon)
            delete(lat)
        end if

        if ($VD .and. hasVoronoi) then
            ; Plot Voronoi cell edges
            system("echo '  Plotting Voronoi cells ...'")
            edgeRes = True
            edgeRes@gsLineThicknessF = 0.5
            edgeRes@gsLineColor = "red"

            do i = 0, numPnt-1
                lon = new(f->numVtx(i)+1, "float")
                lat = new(f->numVtx(i)+1, "float")
                do j = 0, f->numVtx(i)-1
                    lon(j) = lonVtx(i,j)
                    lat(j) = latVtx(i,j)
                end do
                lon(f->numVtx(i)) = lon(0)
                lat(f->numVtx(i)) = lat(0)
                gsn_polyline(wks, map, lon, lat, edgeRes)
                delete(lon)
                delete(lat)
            end do
        end if

        if ($DVT) then
            ; Plot points
            system("echo '  Plotting vertices ...'")
            pntRes = True
            pntRes@gsMarkerIndex = 1
            pntRes@gsMarkerSizeF = 0.005
            pntRes@gsMarkerColor = "red"

            gsn_polymarker(wks, map, lonPnt, latPnt, pntRes)
        end if

        system("echo '  DONE!'")

        frame(wks)
    end
EOF
    if [ -f "$figure.pdf" -a $figtype == "pdf" ]; then
        echo "$figure.pdf generated."
    fi
done
