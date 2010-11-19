#!/bin/bash

echo "Input the data file or wildcard:"
read -p "> " fileWildCard

echo "Plot Delaunay triangulation? (y/n)"
read -p "> " ans
if [ $ans == "y" ]; then
    DT=True
else
    DT=False
fi

echo "Plot Voronoi diagram? (y/n)"
read -p "> " ans
if [ $ans == "y" ]; then
    VD=True
else
    VD=False
fi

echo "Plot vertices? (y/n)"
read -p "> " ans
if [ $ans == "y" ]; then
    DVT=True
else
    DVT=False
fi

mapProj="CE"
echo "Input map projection (NH/SH/CE):"
read -p "[Default: CE] > " ans
if [ -n "$ans" ]; then
    mapProj=$ans
fi

minLat="80"
if [ $mapProj == "NH" ]; then
    echo "Input the southest latitude to view:"
    read -p "[Default: 80] > " ans
    if [ -n "$ans" ]; then
        minLat=$ans
    fi
fi

maxLat="-80"
if [ $mapProj == "SH" ]; then
    echo "Input the northest latitude to view:"
    read -p "[Default: -80] > " ans
    if [ -n "$ans" ]; then
        maxLat=$ans
    fi
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

        numPoint = dimsizes(f->lonPoint)

        lonDVT = f->lonPoint*Rad2Deg
        latDVT = f->latPoint*Rad2Deg

        lonVVT = f->lonVVT*Rad2Deg
        latVVT = f->latVVT*Rad2Deg

        numTriangle = dimsizes(f->triangle(:,0))

        wks = gsn_open_wks("pdf", "$figure")

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

        map = gsn_csm_map(wks, mapRes)

        if (False) then
            ; Plot vertex indices
            textRes = True
            textRes@txFontHeightF = 0.01
            textRes@txFont = "helvetica-bold"

            do i = 0, numPoint-1
                text = gsn_add_text(wks, map, sprinti("%d", i+1), \
                    lonDVT(i)-0.5, latDVT(i)-0.5, textRes)
            end do
        end if

        draw(map)

        if ($DT) then
            ; Plot Delaunay triangle edges
            system("echo Plotting triangle edges ...")
            
            edgeRes = True
            edgeRes@gsLineThicknessF = 0.5
            edgeRes@gsLineColor = "blue"

            lonVtx = new(4, "float")
            latVtx = new(4, "float")
            do i = 0, numTriangle-1
                numVertex = 0
                do j = 0, 2
                    k = f->triangle(i,j)-1
                    if (k .ge. 0) then
                        lonVtx(numVertex) = lonDVT(k)
                        latVtx(numVertex) = latDVT(k)
                        numVertex = numVertex+1
                    end if
                end do
                if (numVertex .gt. 1) then
                    lonVtx(numVertex) = lonVtx(0)
                    latVtx(numVertex) = latVtx(0)
                    gsn_polyline(wks, map, lonVtx(0:numVertex), latVtx(0:numVertex), edgeRes)
                end if
            end do

            delete(lonVtx)
            delete(latVtx)
        end if

        if (False) then
            ; Plot circumcirlces
            system("echo Plotting circumcirlces")

            circRes = True
            circRes@gsLineThicknessF = 0.5
            circRes@gsLineColor = "green"

            arclon = new(50, float)
            arclat = new(50, float)
            do i = 0, numTriangle-1
                nggcog(f->latCirc(i), f->lonCirc(i), f->radius(i), arclat, arclon)
                gsn_polyline(wks, map, arclon, arclat, circRes)
            end do
        end if

        if ($VD) then
            ; Plot Voronoi cell edges
            system("echo Plotting Voronoi cells")

            edgeRes = True
            edgeRes@gsLineThicknessF = 0.5
            edgeRes@gsLineColor = "red"

            do i = 0, numPoint-1
                lonVtx = new(f->numVVT(i)+1, "float")
                latVtx = new(f->numVVT(i)+1, "float")
                do j = 0, f->numVVT(i)-1
                    lonVtx(j) = lonVVT(i,j)
                    latVtx(j) = latVVT(i,j)
                end do
                lonVtx(f->numVVT(i)) = lonVtx(0)
                latVtx(f->numVVT(i)) = latVtx(0)
                gsn_polyline(wks, map, lonVtx, latVtx, edgeRes)
                delete(lonVtx)
                delete(latVtx)
            end do
        end if

        if ($DVT) then
            ; Plot vertices
            system("echo Plotting vertices")

            vertexRes = True
            vertexRes@gsMarkerIndex = 1
            vertexRes@gsMarkerSizeF = 0.005
            vertexRes@gsMarkerColor = "red"

            gsn_polymarker(wks, map, lonDVT, latDVT, vertexRes)
        end if

        frame(wks)

    end
EOF
    if [ -f "$figure.pdf" ]; then
        echo "$figure.pdf generated."
    fi
done
