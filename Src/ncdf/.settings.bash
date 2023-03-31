function var_N {
local var="$1"
case $var in
VAR) printf '%b' 0 ;;
V) printf '%b' 0 ;;
a) printf '%b' 0 ;;
s) printf '%b' 3 ;;
d) printf '%b' 3 ;;
c) printf '%b' 3 ;;
z) printf '%b' 3 ;;
b) printf '%b' 3 ;;
h) printf '%b' 3 ;;
i) printf '%b' 3 ;;
l) printf '%b' 3 ;;
esac
}
