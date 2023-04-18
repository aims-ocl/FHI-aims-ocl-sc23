// Use stdin and stdout as transport method for IPC.

#include "ipc.h"

#define SIGNATURE_STRING "TK_INTERFACE"

int writeSignature()
{
    fprintf(stdout, "%s\n", SIGNATURE_STRING);
    //fprintf(stderr, "writeSignature: %s\n", SIGNATURE_STRING);

    return 0;
}

int readSignature()
{
    char str[13];
    int i_pos;

    i_pos = 0;
    while (1) {
        fscanf(stdin, "%c", &str[i_pos]);
//        fprintf(stderr, "i_pos %i --- %c \n",i_pos,str[i_pos]);

        if (strncmp(str, SIGNATURE_STRING, i_pos + 1)) {
            str[i_pos + 1] = 0;
            i_pos = 0;
            continue;
        }

        i_pos++;

        if (i_pos == strlen(SIGNATURE_STRING)) {
            break;
        }
    }

    return 0;
}

int writeIntegerStdout(const int an_integer)
{
    writeSignature();
    fprintf(stdout, "%i\n", an_integer);
//    fprintf(stderr, "IPC_int (write): %i\n", an_integer);
    fflush (stdout);

    return 1;
}

int readIntegerStdin()
{
    int an_integer;

    readSignature();

    fscanf(stdin, "%i", &an_integer);
//    fprintf(stderr, "IPC_int (read): %i\n", an_integer);
    return an_integer;
}

int writeDoubleArrayStdout(const int n_doubles, const double *double_array)
{
    int i_run = 0;

    for (i_run = 0; i_run < n_doubles; ++i_run) {
        writeSignature();
        fprintf(stdout, "%lf\n", double_array[i_run]);
//		fprintf(stderr, "%i/%i.) writing: %lf\n", i_run, n_doubles, double_array[i_run]);
        fflush (stdout);
    }

    return i_run;
}

int readDoubleArrayStdin(int n_amount, double *doubles)
{

    int i_read;

    for (i_read = 0; i_read < n_amount; ++i_read) {
        readSignature();
        fscanf(stdin, "%lf", &doubles[i_read]);
//        fprintf(stderr, "%i/%i.) IPC_dbl: %lf\n", i_read, n_amount, doubles[i_read]);
    }

    return i_read;
}

void setupStdInOutIPC()
{
    registerWriteInteger(&writeIntegerStdout);
    registerReadInteger(&readIntegerStdin);
    registerWriteDoubleArray(&writeDoubleArrayStdout);
    registerReadDoubleArray(&readDoubleArrayStdin);
}
