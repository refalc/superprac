#include <iostream>
#include "definition.h"

int main(int argc, char *argv[])
{
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if( argc != 3 ) {
        if( rank == 0 ) {
            printf("Wrong number of parameters in command line.\n"
                   "Usage: <ProgName> <Nodes number on (0x) axis> <Nodes number on (0y) axis>\n");
        }

        MPI_Finalize();
        return 1;
    }

    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);

    if( (NY <= 0) || (NX <= 0) ) {
        if (rank == 0)
            printf("The first and the second arguments (mesh numbers) should be positive.\n");

        MPI_Finalize();
        return 2;
    }

    double starttime, endtime;
    starttime = MPI_Wtime();

    Rectangle area;
    area.m_LTPoint.m_dblCrdX = 0.;
    area.m_LTPoint.m_dblCrdY = 0.;
    area.m_RBPoint.m_dblCrdX = 2.;
    area.m_RBPoint.m_dblCrdY = 2.;
    MainFunction(NX, NY, area);

    endtime = MPI_Wtime();
    //printf("Execution time: %f\n", endtime - starttime);

    MPI_Finalize();

    return 0;
}
