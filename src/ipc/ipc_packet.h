// Main purpose is to declare the arrayPacket and all its related functions (creation, destruction,
// querying). The arrayPacket is the quantity sent for communication.

#ifndef IPC_PACKET_H
#define IPC_PACKET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ------------------------------------------------------------------------------------------------
// Declarations
// ------------------------------------------------------------------------------------------------
// Packet type constants. An arrayPacket may carry doubles, integers or strings.
#define PACKET_TYPE_DOUBLE  4001
#define PACKET_TYPE_INTEGER 4002
#define PACKET_TYPE_STRING  4003

// Declare the packet used in the transfer layer.
struct arrayPacket {
    int type_id;
    int i_dimensions;
    int *i_dimension_lengths;
    double *data;
};


// ------------------------------------------------------------------------------------------------
// Creation and destruction.
// ------------------------------------------------------------------------------------------------
// Create an arrayPacket.
struct arrayPacket *createArrayPacket(const int type_id,
                                      const int i_dimensions,
                                      const int *i_dimension_lengths);

// Create an arrayPacket that carries doubles.
struct arrayPacket *createArrayPacketDouble(const int i_dimensions,
                                            const int *i_dimension_lengths,
                                            double *data_double);

// Create an arrayPacket that carries integers.
struct arrayPacket *createArrayPacketInteger(const int i_dimensions,
                                             const int *i_dimension_lengths,
                                             int *data_integer);

// Destroy an arrayPacket.
void freeArrayPacket(struct arrayPacket *packet);


// Create an arrayPacket that carries a string.
struct arrayPacket *createArrayPacketString(const char *data_string);


// ------------------------------------------------------------------------------------------------
// Information inquery.
// ------------------------------------------------------------------------------------------------
// Determines the total array size from the lengths of all its dimensions.
int arraySizeFromDimensionLengths(const int i_dimensions,
                                  const int *i_dimension_lengths);

// Determines if a given arrayPacket carries doubles.
int isArrayPacketDoubles(struct arrayPacket *packet);

// Determines if a given arrayPacket carries integers.
int isArrayPacketIntegers(struct arrayPacket *packet);

// Determines if a given arrayPacket carries a string.
int isArrayPacketString(struct arrayPacket *packet);

// Queries the dimension size of an arrayPacket's data.
int arrayPacketDimensionSize(struct arrayPacket *packet, int i_dimension);

// Determines if two packets are equal.
int arePacketsEqual(const struct arrayPacket *packet1,
                    const struct arrayPacket *packet2);


// ------------------------------------------------------------------------------------------------
// Data copy and conversion.
// ------------------------------------------------------------------------------------------------
// Copy the arrayPacket's data to a double array.
int copyArrayPacketData(const struct arrayPacket *packet, double *data);

// Copy and convert the arrayPacket's data to an int array.
int copyArrayPacketDataToInt(const struct arrayPacket *packet, int *data);

int arrayPacketDataToString(struct arrayPacket *packet,
                            const int max_string_length,
                            char *read_string);


// ------------------------------------------------------------------------------------------------
// Output.
// ------------------------------------------------------------------------------------------------
void printArrayPacket(struct arrayPacket *packet);

#endif // IPC_PACKET_H