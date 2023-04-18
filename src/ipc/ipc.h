// Main purpose is to declare communication commands and functions for communication, i.e. the 
// exchange of strings and arrays with another IPC peer.
#ifndef IPC_H
#define IPC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ipc_packet.h"

// ------------------------------------------------------------------------------------------------
// Function pointer type definitions.
// ------------------------------------------------------------------------------------------------
// Write and read functions.
// A type for a function that can write an integer.
typedef int(*write_integer_function) (const int);
// A type for a function that can read an integer.
typedef int(*read_integer_function)();
// A type for a function that can write a double array.
typedef int(*write_double_array_function)(const int n_amount, const double *double_array);
// A type for a function that can read a double array.
typedef int(*read_double_array_function)(const int n_amount, double *double_array);

// Transaction handler functions.
// A type for a function handling a transaction.
typedef int(*on_transaction)(void);
// A type for a function handling a transaction request.
typedef int(*on_transaction_request)(const char *transaction_name, on_transaction *transaction_handler);

// ------------------------------------------------------------------------------------------------
// Declare transfer layer constants.
// ------------------------------------------------------------------------------------------------
#define SESSION_START 1001
#define SESSION_END   1002

#define TRANSACTION_START      2000
#define TRANSACTION_UNHANDELED 2001
#define TRANSACTION_REJECTED   2002
#define TRANSACTION_ACTIVE     2003

#define WRITE_PACKET_START                  3001
#define WRITE_PACKET_FINISH                 3002
#define READ_PACKET_CONFIRM_START           3003
#define READ_PACKET_CONFIRM_TYPEID          3004
#define READ_PACKET_CONFIRM_DIMENSIONS      3005
#define READ_PACKET_CONFIRM_DIMENSIONLENGTH 3006
#define READ_PACKET_CONFIRM_DATA            3007
#define READ_PACKET_CONFIRM_FINISH          3008

// ------------------------------------------------------------------------------------------------
// Registration of communication carrier. (could be stdout/stdin, sockets, thread safe queues) 
// ------------------------------------------------------------------------------------------------
void registerWriteInteger(write_integer_function function_to_call);
void registerReadInteger(read_integer_function function_to_call);
void registerWriteDoubleArray(write_double_array_function function_to_call);
void registerReadDoubleArray(read_double_array_function function_to_call);


// ------------------------------------------------------------------------------------------------
// Write and read.
// ------------------------------------------------------------------------------------------------
// private internals 
// Writes an integer.
int writeInteger(const int an_integer);
// Reads an integer.
int readInteger();

// Writes a double array.
int writeDoubleArray(const int n_amount, const double *double_array);
// Reads a double array.
int readDoubleArray(const int n_amount, double *double_array);

// Writes an arrayPacket.
int writePacket(struct arrayPacket *packet);
// Reads an arrayPacket.
struct arrayPacket *readPacket();

// For public use. 
// Writes a string.
void writeString(const char *string_to_write);
// Reads a string.
int readString(int max_string_length, char *read_string);

// Writes an array.
void writeArray(const char *array_name, struct arrayPacket *packet);
// Reads an array.
struct arrayPacket *readArray(const char *array_name);


// ------------------------------------------------------------------------------------------------
// Transaction
// ------------------------------------------------------------------------------------------------
// Starts a transaction. The other IPC peer must listen for transactions.
int startTransaction(const char *transaction_name);
// Accept/listen for transaction.
int acceptTransaction(on_transaction_request transaction_registry);


// ------------------------------------------------------------------------------------------------
// Session
// ------------------------------------------------------------------------------------------------
// Ends the session. The other IPC peer must listen for transactions.
void endSession();


// ------------------------------------------------------------------------------------------------
// Private internals
// ------------------------------------------------------------------------------------------------
// Determines if communication is initialized.
int isCommunicationInitialized();

// Checks for expected reply.
void checkReply(const int code_expected,
		        const int reply);

#endif // IPC_H