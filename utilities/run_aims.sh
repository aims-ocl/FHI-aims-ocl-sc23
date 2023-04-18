#!/bin/bash

usage () {
    echo 1>&2 "Usage: $0 [options] BINARY [DIR]

Convenience tool to start FHI-aims in the one-run one-dir paradigma.
Copies control.in and geometry.in into DIR and runs aims on them.

Options:
     -h   --help                 This message
     -d   --dry-run              Dry run
     -c   --control CONTROL      Use this instead of control.in
     -g   --geometry GEOMETRY    Use this instead of geometry.in
     -o   --output OUTFILE       Use this instead of aims.out
     -C   --sed-control SCRIPT   Apply sed-script to control.in.
     -G   --sed-geometry SCRIPT  Apply sed-script to geometry.in.
     -t   --tee                  Tee aims output to stdout
     -n   --mpi-np  NP           Number of CPUs
     -s   --serial               Serial run
     -f   --force                Remove DIR if it exists

Example:
  # Run FHI-aims on k-grids kxkxk with k in 5, 7, 9, 11.
  for k in 5 7 9 11; do
    run_aims.sh -n 2 aims.042811.mpi.x run-$k -C 's/k_grid.*/k_grid $k $k $k/'
  done

Environment:
  AIMS_NEEDS_MPD                 If set (to 'mpd'), make sure it runs.
  AIMS_DEF_NP                    Default for -n.
"
# Undocumented because of limited use without the modcont script:
#    -M   --mod-control ARGS     Call modcont with ARGS on control.in

    exit $1
}
doerror-usage () {
    echo >&2 "$*"
    usage 2
}
doerror () {
    echo >&2 "$*"
    exit 1
}

# ========== Read options

BINARY=
DIR=
CONTROL=control.in
GEOMETRY=geometry.in
OUTPUT=aims.out
SED_CONTROL=
SED_GEOMETRY=
MODCONTS=()
NP="$AIMS_DEF_NP"
NEEDS_MPD="$AIMS_NEEDS_MPD"
TEE=
DRY=
FORCE=
while [ $# -gt 0 ]
  do
  opt=$1
  shift
  case $opt in
      (-h|--help)     usage 0;;
      (-d|--dry-run)  DRY=1;;
      (-c|--control)  CONTROL="$1"; shift;;
      (-g|--geometry) GEOMETRY="$1"; shift;;
      (-o|--output)   OUTPUT="$1"; shift;;
      (-C|--sed-control)   SED_CONTROL+="$1;"; shift;;
      (-G|--sed-geometry)  SED_GEOMETRY+="$1;"; shift;;
      (-M|--mod-control)   MODCONTS+=("$1"); shift;;
      (-n|--mpi-np)   NP="$1"; shift;;
      (-s|--serial)   NP=serial;;
      (-t|--tee)      TEE=1;;
      (-f|--force)    FORCE=1;;
      (-?*)           doerror-usage "Unknown option $opt";;
      (*)
      if [ -z "$BINARY" ]; then
	  BINARY="$opt"
      elif [ -z "$DIR" ]; then
	  DIR="$opt"
      else
	  doerror-usage "Too many arguments"
      fi
  esac
done

if [ -z "$BINARY" ]; then
    doerror-usage "No BINARY specified"
fi
if [ -z "$NP" ]; then
    doerror-usage "Number of processes not specified"
fi
if [ -n "$TEE" ]; then
    TEEOUT=/dev/stdout
else
    TEEOUT=/dev/null
fi


# From now on, immediately exit if something goes wrong:
set -e

# ========== Copy || Edit input files

sed_filter () {
    SCRIPT="$1"; IN="$2"; OUT="$3"
    if [ -n "$SCRIPT" ]; then
	echo "sed $SCRIPT <$IN >$OUT"
	sed "$SCRIPT" <"$IN" >"$OUT"
    else
	cp $IN $OUT
    fi
}
if [ -n "$DIR" ]; then
    if [ -e "$DIR" ]; then
	if [ -n "$FORCE" ]; then
	    rm -rf "$DIR"
	else
	    doerror "Directory '$DIR' already exists."
	fi
    fi
    mkdir -p "$DIR"
    sed_filter "$SED_CONTROL" $CONTROL $DIR/control.in
    sed_filter "$SED_GEOMETRY" $GEOMETRY $DIR/geometry.in
    for ((i=0; i<${#MODCONTS[@]}; i++)); do
	modcont -f $DIR/control.in ${MODCONTS[$i]}
    done
else
    test -e aims.out && doerror "No DIR given, ./aims.out already exists."
    if [ -n "$SED_CONTROL" -o -n "$SED_GEOMETRY" -o "${#MODCONTS[@]}" -gt 0 ]
    then doerror "Cannot edit for in-situ (empty DIR) run."
    fi
    DIR=.
fi

# ========== Run

if [ -n "$DRY" ]; then exit 0; fi

if [ "$NP" = serial ]; then
    echo "nohup $BINARY </dev/null 2>&1 | tee $OUTPUT >$TEEOUT"
    (cd $DIR; set -x; nohup $BINARY </dev/null 2>&1 | tee $OUTPUT >$TEEOUT)
else
    if [ -n "$NEEDS_MPD" ]; then
	# Start mpd if needed
	if ! ps aux |grep -v grep | grep -q "\b$NEEDS_MPD\b"; then
	    echo >&2 "Starting $NEEDS_MPD"
	    nohup $NEEDS_MPD </dev/null >/dev/null 2>/dev/null &
	    sleep 1
	fi
    fi
    echo "nohup mpiexec -n $NP $BINARY </dev/null 2>&1 | tee $OUTPUT >$TEEOUT"
    (cd $DIR;
	nohup mpiexec -n $NP $BINARY </dev/null 2>&1 | tee $OUTPUT >$TEEOUT)
fi
