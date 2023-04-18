#ifndef IPC_STDIO_H
#define IPC_STDIO_H


#include "ipc.h"

int writeSignature();
int readSignature();
int writeIntegerStdout(const int an_integer);
int readIntegerStdin();
int writeDoubleArrayStdout(int n_amount, double *doubles);
int readDoubleArrayStdin(int n_amount, double *doubles);
void setupStdInOutIPC();

#endif