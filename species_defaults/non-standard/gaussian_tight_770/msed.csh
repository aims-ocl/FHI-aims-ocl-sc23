#!/bin/csh

foreach x ($argv)
    sed s/'basis_dep_cutoff    1e-4'/'basis_dep_cutoff    0e-0'/g $x > test.tt
    mv test.tt $x
end
