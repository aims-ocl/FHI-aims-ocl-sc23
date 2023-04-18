#!/usr/bin/awk -f

BEGIN {
    print "# aims-extract_input_files.awk (j.meyer@lic.leidenuniv.nl)"
}

function print_lines_between_separator(filename) {

    separator = "-----------------------------------------------------------------------"

    do {
        getline
    } while ($1 != separator)

    getline ; getline

    while ($1 != separator) {
#        print substr($0,3)
        printf substr($0,3) "\n" >> filename
        getline
    }
}

/Parsing control.in/ {

    control_in_file = FILENAME "-control.in"
    print_lines_between_separator(control_in_file)
    printf "|-> %s written \n", control_in_file
}

/Parsing geometry.in/ {

    geometry_in_file = FILENAME "-geometry.in"
    print_lines_between_separator(geometry_in_file)
    printf "|-> %s written \n", geometry_in_file
}
