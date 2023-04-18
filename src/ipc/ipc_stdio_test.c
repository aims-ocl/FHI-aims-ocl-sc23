// Unittest for ipc_stdio

#include "ipc_stdio.h"
#include <string.h>
#include <assert.h>

enum ROLE {
    MASTER, SLAVE, UNSET
};

double sumDoubles(int n_doubles, double *data)
{
    int i_run;
    double result = 0.0;

    for (i_run = 0; i_run < n_doubles; ++i_run) {
        result += data[i_run];
    }

    return result;
}

void testCommunicateIntegerMaster()
{
    fprintf(stderr, "Test communicate integer (MASTER)\n");

    int an_integer = 10, reply;

    writeIntegerStdout(an_integer);
    reply = readIntegerStdin();

    assert(reply == an_integer);

    fprintf(stderr, "MASTER done\n");
}

void testCommunicateIntegerSlave()
{
    fprintf(stderr, "Test communicate integer (SLAVE)\n");

    int an_integer;

    an_integer = readIntegerStdin();
    writeIntegerStdout(an_integer);

    fprintf(stderr, "SLAVE done\n");
}

void testCommunicateDoubleArrayMaster()
{
    fprintf(stderr, "Test communicate double array (MASTER)\n");

    const int n_doubles = 5000;
    int i_run;
    double data_doubles[n_doubles], reply_sum, result_diff;

    for (i_run = 0; i_run < n_doubles; ++i_run) {
        data_doubles[i_run] = i_run * 1.1;
    }

    writeIntegerStdout(n_doubles);
    writeDoubleArrayStdout(n_doubles, data_doubles);
    readDoubleArrayStdin(1, &reply_sum);

    result_diff = sumDoubles(n_doubles, data_doubles) - reply_sum;

    fprintf(stderr,
            "reply sum: %lf --- Result diff %lf\n",
            reply_sum,
            result_diff);

    assert(abs(result_diff) < 1.E-8);

    fprintf(stderr, "MASTER done\n");
}

void testCommunicateDoubleArraySlave()
{
    fprintf(stderr, "Test communicate double array (SLAVE)\n");

    int n_doubles, i_read;
    double *data_double, sum;

    n_doubles = readIntegerStdin();

    data_double = malloc(n_doubles * sizeof(double));

    i_read = readDoubleArrayStdin(n_doubles, data_double);

    assert(i_read == n_doubles);

    sum = sumDoubles(n_doubles, data_double);
    free(data_double);

    writeDoubleArrayStdout(1, &sum);

    fprintf(stderr, "SLAVE done\n");
}

int main(int argc, char *argv[])
{
    enum ROLE role = UNSET;

    fprintf(stderr, "Starting test.\n");

    if (argc > 1) {
        if (!strcmp(argv[1], "master")) {
            role = MASTER;
        }
        else if (!strcmp(argv[1], "slave")) {
            role = SLAVE;
        }
    }

    if (role == MASTER) {
        testCommunicateIntegerMaster();
        testCommunicateDoubleArrayMaster();
    }

    if (role == SLAVE) {
        testCommunicateIntegerSlave();
        testCommunicateDoubleArraySlave();
    }

    exit(0);
}
