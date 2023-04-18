// ArrayPackets are used to exchange data between two IPC peers. 

#include "ipc_packet.h"

// Determine total array size from its dimension sizes.
int arraySizeFromDimensionLengths(const int i_dimensions,
                                  const int *i_dimension_lengths)
{
    int i_run, array_size;

    array_size = 1;
    for (i_run = 0; i_run < i_dimensions; ++i_run) {
        array_size *= i_dimension_lengths[i_run];
    }

    return array_size;
}

// Create an arrayPacket.
struct arrayPacket *createArrayPacket(const int type_id,
                                      const int i_dimensions,
                                      const int *i_dimension_lengths)
{
    struct arrayPacket *packet;
    int array_size;
    int i_run;

    packet = malloc(sizeof(struct arrayPacket));
    if(packet == NULL){
        printf("Failed to alloc arrayPacket structure.\n");
        return NULL;
    }

    packet->type_id = type_id;
    packet->i_dimensions = i_dimensions;
    packet->i_dimension_lengths = malloc(i_dimensions * sizeof(int));
    if(packet->i_dimension_lengths == NULL){
        printf("Failed to alloc dimension length.\n");
        return NULL;
    }
    

    for (i_run = 0; i_run < i_dimensions; ++i_run) {
        packet->i_dimension_lengths[i_run] = i_dimension_lengths[i_run];
    }

    array_size = arraySizeFromDimensionLengths(i_dimensions,
                                               i_dimension_lengths);

    packet->data = malloc(array_size * sizeof(double));
    if(packet->data == NULL){
        printf("Failed to alloc packet data.\n");
        return NULL;
    }
    return packet;
}

// Create arrayPacket from doubles.
struct arrayPacket *createArrayPacketDouble(const int i_dimensions,
                                            const int *i_dimension_lengths,
                                            double *data_double)
{
    struct arrayPacket *packet;
    const int type_id = PACKET_TYPE_DOUBLE;
    const int array_size = arraySizeFromDimensionLengths(i_dimensions,
                                                         i_dimension_lengths);

    packet = createArrayPacket(type_id, i_dimensions, i_dimension_lengths);
    memcpy(packet->data, data_double, array_size * sizeof(double));
    return packet;
}

// Create arrayPacket from integers.
struct arrayPacket *createArrayPacketInteger(const int i_dimensions,
                                             const int *i_dimension_lengths,
                                             int *data_integer)
{
    struct arrayPacket *packet;
    const int type_id = PACKET_TYPE_INTEGER;
    const int array_size = arraySizeFromDimensionLengths(i_dimensions,
                                                         i_dimension_lengths);

    int i_run;

    packet = createArrayPacket(type_id, i_dimensions, i_dimension_lengths);

    for (i_run = 0; i_run < array_size; ++i_run) {
        packet->data[i_run] = data_integer[i_run];
    }
    return packet;
}

// Create arrayPacket from a string.
struct arrayPacket *createArrayPacketString(const char *data_string)
{
    struct arrayPacket *packet;
    const int i_dimensions = 1;
    int i_dimension_lengths[1];
    const int type_id = PACKET_TYPE_STRING;
    int i_index;

    i_dimension_lengths[0] = strlen(data_string);
    packet = createArrayPacket(type_id, i_dimensions, i_dimension_lengths);

    for (i_index = 0; i_index < i_dimension_lengths[0]; ++i_index) {
        packet->data[i_index] = data_string[i_index];
    }
    return packet;
}

// Determines if the arrayPacket carries double data.
int isArrayPacketDoubles(struct arrayPacket *packet){
    if(packet->type_id == PACKET_TYPE_DOUBLE) {
        return 1;
    }
    
    return 0;
}

// Determines if the arrayPacket carries integer data.
int isArrayPacketIntegers(struct arrayPacket *packet){
    if(packet->type_id == PACKET_TYPE_INTEGER) {
        return 1;
    }

    return 0;
}

// Determines if the arrayPacket carries string data.
int isArrayPacketString(struct arrayPacket *packet){
    if(packet->type_id == PACKET_TYPE_STRING) {
        return 1;
    }

    return 0;
}

// Determines the dimension size of a dimenion of the packet's data.
int arrayPacketDimensionSize(struct arrayPacket *packet, int i_dimension){
    if(i_dimension > packet->i_dimensions){
        return 0;
    }
    
    return packet->i_dimension_lengths[i_dimension-1];
}

// Copies the packet's data to a double array.
int copyArrayPacketData(const struct arrayPacket *packet, double *data){
    int array_size = arraySizeFromDimensionLengths(packet->i_dimensions,
                                                   packet->i_dimension_lengths);
    
    memcpy(data, packet->data, array_size * sizeof(double));

    return array_size;
}

// Copies and converts the packet's data to an integer array.
int copyArrayPacketDataToInt(const struct arrayPacket *packet, int *data){
    int array_size = arraySizeFromDimensionLengths(packet->i_dimensions,
                                                   packet->i_dimension_lengths);
    int i_run;
    
    for(i_run = 0; i_run < array_size; ++i_run){
            data[i_run] = packet->data[i_run]; 
    }
    
    return array_size;
}

// Destroies an arrayPacket.
void freeArrayPacket(struct arrayPacket *packet)
{
    free(packet->i_dimension_lengths);
    free(packet->data);
    free(packet);
}

// Print an arrayPacket.
void printArrayPacket(struct arrayPacket *packet)
{
    int i_run, array_size;
    printf("Printing array packet\n");
    printf("TypeID: %i\n", packet->type_id);
    printf("Dimensions: %i\n", packet->i_dimensions);
    printf("Dimensions lengths: \n");
    for (i_run = 0; i_run < packet->i_dimensions; ++i_run) {
        printf("%i\n", packet->i_dimension_lengths[i_run]);
    }

    array_size = arraySizeFromDimensionLengths(packet->i_dimensions,
                                               packet->i_dimension_lengths);

    printf("Data: \n");
    for (i_run = 0; i_run < array_size; ++i_run) {
        printf("%lf\n", packet->data[i_run]);
    }
}

// Determines if two arrayPackets are equal.
int arePacketsEqual(const struct arrayPacket *packet1,
                    const struct arrayPacket *packet2)
{

    int i_dimension, i_element, array_size;

    if (packet1->type_id != packet2->type_id) {
        return 0;
    }

    if (packet1->i_dimensions != packet2->i_dimensions) {
        return 0;
    }

    for (i_dimension = 0; i_dimension < packet1->i_dimensions; ++i_dimension) {
        if (packet1->i_dimension_lengths[i_dimension]
                != packet2->i_dimension_lengths[i_dimension]) {
            return 0;
        }
    }

    array_size = arraySizeFromDimensionLengths(packet1->i_dimensions,
                                               packet1->i_dimension_lengths);

    for (i_element = 0; i_element < array_size; ++i_element) {
        if (abs(packet1->data[i_element] - packet2->data[i_element]) > 1.0E-8) {
            return 0;
        }
    }

    return 1;
}

// Converts the packet's data to a string.
int arrayPacketDataToString(struct arrayPacket *packet,
                            const int max_string_length,
                            char *read_string)
{

    const int array_size =
            arraySizeFromDimensionLengths(packet->i_dimensions,
                                          packet->i_dimension_lengths);

    int i_run, ascii_code;

    for (i_run = 0; i_run < array_size && i_run < max_string_length-1; ++i_run) {
        ascii_code = packet->data[i_run];
        read_string[i_run] = ascii_code;
    }
    read_string[i_run] = 0;

    return i_run;
}
