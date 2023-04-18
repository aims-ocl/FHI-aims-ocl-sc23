// Allows to send strings and arrays between two IPC peers.

#include "ipc.h"

// Define the low level transfer layer procedure pointer.
// Only these functions need concrete implementatoins for the IPC to work.
// A concrete implementation could use sockets, stdin/stdout, thread safe queues, ...
// These pointers are private and can be set through the corresponding register functions.
write_integer_function write_integer_implementation = NULL;
read_integer_function read_integer_implementation = NULL;
write_double_array_function write_double_array_implementation = NULL;
read_double_array_function read_double_array_implementation = NULL;

// Registers the function to write integers.
void registerWriteInteger(write_integer_function function_to_call)
{
    write_integer_implementation = function_to_call;
}

// Registers the function to read integers.
void registerReadInteger(read_integer_function function_to_call)
{
    read_integer_implementation = function_to_call;
}

// Registers the function to write a double array.
void registerWriteDoubleArray(write_double_array_function function_to_call)
{
    write_double_array_implementation = function_to_call;
}

// Registers the function to read a double array.
void registerReadDoubleArray(read_double_array_function function_to_call)
{
    read_double_array_implementation = function_to_call;
}

// Write an integer to the other IPC peer.
int writeInteger(const int an_integer)
{
    return write_integer_implementation(an_integer);
}

// Read an integer from the other IPC peer.
int readInteger()
{
    return read_integer_implementation();
}

// Write a double array to the other IPC peer.
int writeDoubleArray(const int n_amount, const double *double_array)
{
    return write_double_array_implementation(n_amount, double_array);
}

// Read a double array from the other IPC peer.
int readDoubleArray(const int n_amount, double *double_array)
{
    return read_double_array_implementation(n_amount, double_array);
}

// Write a string to the other IPC peer.
void writeString(const char *string_to_write)
{

    struct arrayPacket *packet;

    packet = createArrayPacketString(string_to_write);
    writePacket(packet);
    freeArrayPacket(packet);
}

// Reads a string from the other IPC peer.
int readString(int max_string_length, char *read_string)
{
    struct arrayPacket *packet;

    int number_read_chars;

    packet = readPacket();
    number_read_chars = arrayPacketDataToString(packet,
                                                max_string_length,
                                                read_string);
    freeArrayPacket(packet);

    return number_read_chars;
}

// Accecpt/listen for a transaction.
int acceptTransaction(on_transaction_request transaction_registry)
{
    int command;
    char transaction_name[1024];
    struct arrayPacket *packet;
    on_transaction *transaction_handler;

    while (1) {
//        printf("About to read command %i\n", command);
        command = readInteger();
//        printf("Read command: %i\n", command);

        if (command == TRANSACTION_START) {
            packet = readArray("START_TRANSACTION");
            arrayPacketDataToString(packet,
                                    sizeof(transaction_name),
                                    transaction_name);

//            printf("Invoking transaction registry: %s\n", transaction_name);

            // Ask for transaction handler of this transaction.
            transaction_registry(transaction_name, transaction_handler);

            // If transaction is handled execute the handler. Otherwise inform peer.
            if (transaction_handler != NULL) {
//            printf( "Invoking transaction handler: %s\n", transaction_name);
                writeInteger (TRANSACTION_ACTIVE);
                (*transaction_handler)();
            }
            else {
                writeInteger (TRANSACTION_UNHANDELED);
            }
        }
        else if (command == SESSION_END) {
//            printf("Ending session\n");
            break;
        }
        else {
            printf("Received unexpected command code %i\n. Abort.", command);
            exit(-1);
        }
    }
    return 0;
}

// Determines if communication is initialized.
int isCommunicationInitialized(){
    // Write integer
    if(write_integer_implementation == NULL){
        return 0;
    }
    // Read integer
    if(read_integer_implementation == NULL){
        return 0;
    }
    
    // Write double array
    if(write_double_array_implementation == NULL){
        return 0;
    }
    
    // Read double array
    if(read_double_array_implementation == NULL){
        return 0;
    }

    // Every essential function is set.
    return 1;
}

// Starts a transaction.
int startTransaction(const char *transaction_name)
{

    int reply;

    // Determine if communication is initialized.
    if (isCommunicationInitialized() == 0) {
        return 0;
    }
    
    // Command to start transaction. The name of the transaction.
    writeInteger (TRANSACTION_START);
    writeString("START_TRANSACTION");
    writeString(transaction_name);

    // Read transaction state.
    reply = readInteger();

    // Unhandeled of rejected.
    if (reply == TRANSACTION_UNHANDELED || reply == TRANSACTION_REJECTED) {
        return 0;
    }

    // Transaction active.
    if (reply == TRANSACTION_ACTIVE) {
        return 1;
    }

    printf("Unexpected and unknown transaction state.");
    exit(-1);
}

// Ends a session. IPC peer must listen for transaction.
void endSession()
{
    if (write_integer_implementation == NULL) {
        return;
    }

    writeInteger (SESSION_END);
}

// Checks that a code/reply equals an expected one.
void expectCode(const int code_expected, const int reply)
{
    if (code_expected != reply) {
        printf("Got code %i but expected %i. Program abort.\n",
               code_expected,
               reply);
        exit(-1);
    }
}

// Writes an arrayPacket.
int writePacket(struct arrayPacket *packet)
{
    int reply, array_size, i_dimension;

    writeInteger (WRITE_PACKET_START);
    reply = readInteger();
    expectCode(READ_PACKET_CONFIRM_START, reply);

    writeInteger(packet->type_id);
    reply = readInteger();
    expectCode(READ_PACKET_CONFIRM_TYPEID, reply);

    writeInteger(packet->i_dimensions);
    reply = readInteger();
    expectCode(READ_PACKET_CONFIRM_DIMENSIONS, reply);

    for (i_dimension = 0; i_dimension < packet->i_dimensions; ++i_dimension) {
        writeInteger(packet->i_dimension_lengths[i_dimension]);
        reply = readInteger();
        expectCode(READ_PACKET_CONFIRM_DIMENSIONLENGTH, reply);
    }

    array_size = arraySizeFromDimensionLengths(packet->i_dimensions,
                                               packet->i_dimension_lengths);
    writeDoubleArray(array_size, packet->data);

    writeInteger (WRITE_PACKET_FINISH);
    reply = readInteger();
    expectCode(READ_PACKET_CONFIRM_FINISH, reply);

    return reply;
}

// Reads an arrayPacket.
struct arrayPacket *readPacket()
{
    int operation, array_size, i_dimension, i_dimensions, type_id;
    int *i_dimension_lengths;
    struct arrayPacket *packet;

    operation = readInteger();
    expectCode(WRITE_PACKET_START, operation);
    writeInteger (READ_PACKET_CONFIRM_START);

    type_id = readInteger();
    writeInteger (READ_PACKET_CONFIRM_TYPEID);

    i_dimensions = readInteger();
    writeInteger (READ_PACKET_CONFIRM_DIMENSIONS);

    i_dimension_lengths = malloc(i_dimensions * sizeof(int));

    for (i_dimension = 0; i_dimension < i_dimensions; ++i_dimension) {
        i_dimension_lengths[i_dimension] = readInteger();
        writeInteger (READ_PACKET_CONFIRM_DIMENSIONLENGTH);
    }

    packet = createArrayPacket(type_id, i_dimensions, i_dimension_lengths);

    array_size = arraySizeFromDimensionLengths(packet->i_dimensions,
                                               packet->i_dimension_lengths);

    readDoubleArray(array_size, packet->data);

    operation = readInteger();
    expectCode(WRITE_PACKET_FINISH, operation);
    writeInteger (READ_PACKET_CONFIRM_FINISH);
    free(i_dimension_lengths);

    return packet;
}

// Writes an array.
void writeArray(const char *array_name, struct arrayPacket *packet)
{
    writeString(array_name);
    writePacket(packet);
}

// Reads an array.
struct arrayPacket *readArray(const char *array_name)
{

    char array_name_read[1024];
    int number_read_chars;
    struct arrayPacket *packet;

    number_read_chars = readString(sizeof(array_name_read), array_name_read);

    if (strncmp(array_name_read, array_name, number_read_chars)!=0) {
        printf("Number chars unequal %i %c\n", number_read_chars, array_name_read[number_read_chars]);
        return NULL;
    }

    packet = readPacket();

    return packet;
}

