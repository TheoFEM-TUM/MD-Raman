#!/usr/bin/gawk -f

BEGINFILE { start = 0 }

/DIELECTRIC TENSOR \(including/ {
	start = start ? start : FNR
}

(start) && (FNR - start == 2) { xx = $1; xy = $2; xz = $3 }
(start) && (FNR - start == 3) { yx = $1; yy = $2; yz = $3 }
(start) && (FNR - start == 4) { zx = $1; zy = $2; zz = $3 }

ENDFILE {
	printf "  %+12.6f  %+12.6f  %+12.6f", xx, yy, zz
	printf "  %+12.6f  %+12.6f  %+12.6f", xy, yz, zx
	printf "\n"
}
