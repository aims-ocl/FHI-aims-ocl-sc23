// Unittests for ipc_packet

#include <assert.h>
#include "ipc_packet.h"

void testArraySizeFromDimensionLengths()
{
    const int dimensions = 2;
    const int dimension_lengths[2] = { 10, 5 };

    int array_size;

    printf("Tests arraySizeFromDimensionLengths\n");

    array_size = arraySizeFromDimensionLengths(dimensions, dimension_lengths);

    assert(array_size == 50);
}

void testCreateArrayPacket()
{
    const int type_id = PACKET_TYPE_DOUBLE;
    const int dimensions = 2;
    const int dimension_lengths[2] = { 3, 10 };

    struct arrayPacket *packet;

    printf("Tests createArrayPacket\n");

    packet = createArrayPacket(type_id, dimensions, dimension_lengths);

    assert(packet->type_id == type_id);
    assert(packet->i_dimensions == dimensions);
    assert(packet->i_dimension_lengths[0] == dimension_lengths[0]);
    assert(packet->i_dimension_lengths[1] == dimension_lengths[1]);

    freeArrayPacket(packet);
}

void testCreateArrayPacketDouble()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    struct arrayPacket *packet;

    printf("Tests createArrayPacketDouble\n");

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);

    assert(packet->type_id == PACKET_TYPE_DOUBLE);
    assert(packet->i_dimensions == dimensions);
    assert(packet->i_dimension_lengths[0] == dimension_lengths[0]);
    assert(packet->i_dimension_lengths[1] == dimension_lengths[1]);

    assert(memcmp(packet->data, data_double, 30 * sizeof(double)) == 0);

    freeArrayPacket(packet);
}

void testCreateArrayPacketInteger()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 2, 10 };
    int data_integer[2 * 10];
    int i_run, i_data;
    struct arrayPacket *packet;

    printf("Tests createArrayPacketInteger\n");

    packet = createArrayPacketInteger(dimensions,
                                      dimension_lengths,
                                      data_integer);

    assert(packet->type_id == PACKET_TYPE_INTEGER);
    assert(packet->i_dimensions == dimensions);
    assert(packet->i_dimension_lengths[0] == dimension_lengths[0]);
    assert(packet->i_dimension_lengths[1] == dimension_lengths[1]);

    for (i_run = 0; i_run < 20; ++i_run) {
        i_data = packet->data[i_run];
        assert(i_data == data_integer[i_run]);
    }
    freeArrayPacket(packet);
}

void testCreateArrayPacketString()
{
    const char data_string[] = "A test string";

    struct arrayPacket *packet;

    printf("Tests createArrayString\n");

    packet = createArrayPacketString(data_string);

    assert(1 == 1);

    freeArrayPacket(packet);
}

void testIsArrayPacketDoubles()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    int data_integer[3 * 10];
    struct arrayPacket *packet;

    printf("Tests isArrayPacketDoubles\n");

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);
    assert(isArrayPacketDoubles(packet) == 1);
    freeArrayPacket(packet);

    packet = createArrayPacketInteger(dimensions,
                                      dimension_lengths,
                                      data_integer);
    assert(isArrayPacketDoubles(packet) == 0);
    freeArrayPacket(packet);

    packet = createArrayPacketString("Any string");
    assert(isArrayPacketDoubles(packet) == 0);
    freeArrayPacket(packet);
}

void testIsArrayPacketIntegers()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    int data_integer[3 * 10];
    struct arrayPacket *packet;

    printf("Tests isArrayPacketIntegers\n");

    packet = createArrayPacketInteger(dimensions,
                                      dimension_lengths,
                                      data_integer);
    assert(isArrayPacketIntegers(packet) == 1);
    freeArrayPacket(packet);

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);
    assert(isArrayPacketIntegers(packet) == 0);
    freeArrayPacket(packet);

    packet = createArrayPacketString("Any string");
    assert(isArrayPacketIntegers(packet) == 0);
    freeArrayPacket(packet);
}

void testIsArrayPacketString()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    int data_integer[3 * 10];
    struct arrayPacket *packet;

    printf("Tests isArrayPacketString\n");

    packet = createArrayPacketString("Any string");
    assert(isArrayPacketString(packet) == 1);
    freeArrayPacket(packet);

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);
    assert(isArrayPacketString(packet) == 0);
    freeArrayPacket(packet);

    packet = createArrayPacketInteger(dimensions,
                                      dimension_lengths,
                                      data_integer);
    assert(isArrayPacketString(packet) == 0);
    freeArrayPacket(packet);
}

void testArrayPacketDimensionSize()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    struct arrayPacket *packet;

    printf("Tests arrayPacketDimensionSize\n");

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);
    
    assert(arrayPacketDimensionSize(packet,1) == 3);
    assert(arrayPacketDimensionSize(packet,2) == 10);
    assert(arrayPacketDimensionSize(packet,3) == 0);
    
    freeArrayPacket(packet);
}

void testFreeArrayPacket()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    struct arrayPacket *packet;

    printf("Tests freeArrayPacket\n");

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);

    freeArrayPacket(packet);
}

void testCopyArrayPacketData()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10], data_copied[3 * 10];
    struct arrayPacket *packet;

    int i_run, i_copied;
    
    printf("Tests copyArrayPacketData\n");

    for(i_run = 0; i_run < 30; ++i_run){
        data_double[i_run] = i_run * 1.23;
    }
    
    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);

    i_copied = copyArrayPacketData(packet, data_copied);

    assert(i_copied == 30);
    
    for(i_run = 0; i_run < 30; ++i_run){
        assert(abs(data_double[i_run] - data_copied[i_run]) < 1.e-8);
    }
    
    
    freeArrayPacket(packet);
}

void testCopyArrayPacketDataToInt()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    int data_integer[3 * 10], data_copied[3 * 10];
    struct arrayPacket *packet;

    int i_run, i_copied;
    
    printf("Tests copyArrayPacketDataToInt\n");

    for(i_run = 0; i_run < 30; ++i_run){
        data_integer[i_run] = i_run;
    }
    
    packet = createArrayPacketInteger(dimensions,
                                      dimension_lengths,
                                      data_integer);

    i_copied = copyArrayPacketDataToInt(packet, data_copied);

    assert(i_copied == 30);
    
    for(i_run = 0; i_run < 30; ++i_run){
        assert(data_integer[i_run] == data_copied[i_run]);
    }
    
    freeArrayPacket(packet);
}

void testPrintArrayPacket()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    struct arrayPacket *packet;

    printf("Tests printArrayPacket\n");

    packet = createArrayPacketDouble(dimensions,
                                     dimension_lengths,
                                     data_double);

    printArrayPacket(packet);
    freeArrayPacket(packet);
}

void testArePacketsEqual()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];

    struct arrayPacket *packet1;
    struct arrayPacket *packet2;

    printf("Tests arePacketsEqual\n");

    packet1 = createArrayPacketDouble(dimensions,
                                      dimension_lengths,
                                      data_double);

    packet2 = createArrayPacketDouble(dimensions,
                                      dimension_lengths,
                                      data_double);

    assert(arePacketsEqual(packet1, packet1) == 1);
    assert(arePacketsEqual(packet1, packet2) == 1);

    packet2->data[1] -= 1.0;
    assert(arePacketsEqual(packet1, packet2) == 0);
    freeArrayPacket(packet2);

    packet2 = createArrayPacketDouble(dimensions,
                                      dimension_lengths,
                                      data_double);

    packet2->i_dimensions = 1;
    assert(arePacketsEqual(packet1, packet2) == 0);

    freeArrayPacket(packet2);

    packet2 = createArrayPacketDouble(dimensions,
                                      dimension_lengths,
                                      data_double);

    packet2->i_dimension_lengths[1] = 5;
    assert(arePacketsEqual(packet1, packet2) == 0);

    freeArrayPacket(packet1);
    freeArrayPacket(packet2);
}

void testArrayPacketDataToString()
{

    const char data_string[] = "A test string";

    struct arrayPacket *packet;

    char str_from_packet[80];
    int number_chars_read;

    printf("Tests arrayPacketDataToString\n");

    packet = createArrayPacketString(data_string);

    arrayPacketDataToString(packet, sizeof(str_from_packet), str_from_packet);

    assert(!strcmp(data_string, str_from_packet));

    number_chars_read = arrayPacketDataToString(packet, 5, str_from_packet);

    assert(number_chars_read == 5);
    assert(!strncmp(data_string, str_from_packet, number_chars_read));

    freeArrayPacket(packet);
}

int main(void)
{

    printf("Starting test.\n");

    testArraySizeFromDimensionLengths();
    testCreateArrayPacket();
    testCreateArrayPacketDouble();
    testCreateArrayPacketInteger();
    testCreateArrayPacketString();
    testIsArrayPacketDoubles();
    testIsArrayPacketIntegers();
    testIsArrayPacketString();
    testArrayPacketDimensionSize();
    testCopyArrayPacketData();
    testCopyArrayPacketDataToInt();
    testFreeArrayPacket();
    testPrintArrayPacket();
    testArePacketsEqual();
    testArrayPacketDataToString();
    

    printf(".all tests done\n");
    return 0;
}
