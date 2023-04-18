// Unittest for IPC.

#include "ipc.h"
#include "ipc_stdio.h"
#include <assert.h>

enum ROLE {
    MASTER, SLAVE, UNSET
};

struct arrayPacket *fillTestPacket()
{
    const int dimensions = 2;
    const int dimension_lengths[] = { 3, 10 };
    double data_double[3 * 10];
    struct arrayPacket *test_packet;
    int i_run;

    for (i_run = 0; i_run < (sizeof(data_double) / sizeof(double)); ++i_run) {
        data_double[i_run] = i_run * 1.23;
    }

    test_packet = createArrayPacketDouble(dimensions,
                                          dimension_lengths,
                                          data_double);

    return test_packet;
}

void testCommunicatePacketMaster()
{
    struct arrayPacket *packet, *read_packet;

    fprintf(stderr, "Test communicate packet (MASTER)\n");

    packet = fillTestPacket();
    writePacket(packet);
    read_packet = readPacket();

    assert(arePacketsEqual(packet, read_packet) == 1);

    freeArrayPacket(packet);
    freeArrayPacket(read_packet);

    fprintf(stderr, "MASTER done\n");
}

void testCommunicatePacketSlave()
{
    struct arrayPacket *read_packet;

    fprintf(stderr, "Test communicate packet (SLAVE)\n");

    read_packet = readPacket();

    writePacket(read_packet);

    freeArrayPacket(read_packet);

    fprintf(stderr, "SLAVE done\n");
}

struct arrayPacket *createTransferTestPacket()
{

    int type_id, i_dimensions, i_dimension_lengths[2], i_element;
    struct arrayPacket *packet;

    type_id = 1;
    i_dimensions = 2;
    i_dimension_lengths[0] = 3000;
    i_dimension_lengths[1] = 50000;

    i_dimension_lengths[1] = 50;

    packet = createArrayPacket(type_id, i_dimensions, i_dimension_lengths);

    for (i_element = 0;
            i_element < i_dimension_lengths[0] * i_dimension_lengths[1];
            ++i_element) {
        packet->data[i_element] = i_element;
    }

    return packet;
}

struct arrayPacket *createTransferTestPacketInteger()
{
    int i_dimensions, i_dimension_lengths[1], i_element, i_data[2000];
    struct arrayPacket *packet;

    i_dimensions = 2;
    i_dimension_lengths[0] = 100;
    i_dimension_lengths[1] = 20;

    for (i_element = 0; i_element < sizeof(2000) / sizeof(int); ++i_element) {
        i_data[i_element] = -i_element;
    }

    packet = createArrayPacketInteger(i_dimensions,
                                      i_dimension_lengths,
                                      i_data);

    return packet;
}

void minimalTest()
{
    double *test_array;
    int *test_array_integer;

    int i_index;

    struct arrayPacket *packet;

    printf("FHI-aims: MINIMAL test starts.\n");
    printf("FHI-aims: Start test transaction (READ).\n");

    if (startTransaction("TRANSFER_TEST_READ") == 1) {
        printf("FHI-aims: Create test packet.\n");
        packet = createTransferTestPacket();
        printf("FHI-aims: Sending packet.\n");
        writePacket(packet);
        freeArrayPacket(packet);
        printf("FHI-aims: Sent packet.\n");
    }

    printf("FHI-aims: Start test transaction (READ).\n");
    if (startTransaction("TRANSFER_TEST_READ") == 1) {
        printf("FHI-aims: Create test packet. (integer)\n");
        packet = createTransferTestPacketInteger();
        printf("FHI-aims: Sending packet.\n");
        writePacket(packet);
        freeArrayPacket(packet);
        printf("FHI-aims: Sent packet.\n");
    }

    printf("FHI-aims: Start test transaction (READ).\n");
    if (startTransaction("TRANSFER_TEST_READ") == 1) {
        printf("FHI-aims: Create test packet. (string)\n");
        packet =
                createArrayPacketString("aBcDeFgHiJ DAS ist EIN test STRING ABCbca");
        printf("FHI-aims: Sending packet.\n");
        writePacket(packet);
        freeArrayPacket(packet);
        printf("FHI-aims: Sent packet.\n");
    }

    printf("FHI-aims: Start test transaction (WRITE).\n");
    if (startTransaction("TRANSFER_TEST_WRITE") == 1) {
        printf("FHI-aims: Reading packet.\n");
        packet = readPacket();
        printf("FHI-aims: Read packet.\n");
        printArrayPacket(packet);
    }

    printf("FHI-aims: Sent the recieved packet back.\n");
    if (startTransaction("TRANSFER_TEST_READ") == 1) {
        writePacket(packet);
    }

    if (startTransaction("TRANSFER_TEST_READ_ARRAY") == 1) {
        packet = createArrayPacketString("a string transfered");
        writeArray("TEST_ARRAY_STRING", packet);
        freeArrayPacket(packet);
    }

    endSession();
}

// To execute tests: (assuming test compiled to ./ipc_test)
// >>> mkfifo fifo0 fifo1
// >>> ( ./ipc_test master > fifo0 < fifo1 ) &  ( ./ipc_test slave > fifo1 < fifo0 ) & ( exec 30<fifo0 31<fifo1 )
int main(int argc, char *argv[])
{
    enum ROLE role = UNSET;

    fprintf(stderr, "Starting test.\n");

    setupStdInOutIPC();

    if (argc > 1) {
        if (!strcmp(argv[1], "master")) {
            role = MASTER;
        }
        else if (!strcmp(argv[1], "slave")) {
            role = SLAVE;
        }
    }

    if (role == MASTER) {
        testCommunicatePacketMaster();
    }

    if (role == SLAVE) {
        testCommunicatePacketSlave();
    }

    exit(0);
}
