#include "definition.h"
#include <math.h>
#include <limits>
#include <iostream>
#include <algorithm>

//----MATH FUNCTIONS----
double Boundary(double x, double y)
{
    return 2. / (1. + x*x + y*y);
}

double RightSide(double x, double y)
{
    double top = 8. * (1. - x*x - y*y);
    double bot = (1. + x*x + y*y) * (1. + x*x + y*y) * (1. + x*x + y*y);

    return top / bot;
}

double LaplasOperator(const CMatrix &matrix, const CMesh &mesh, int x, int y)
{
    const double ldx = (matrix(x, y) - matrix(x - 1, y)) / mesh.GetStepX();
    const double rdx = (matrix(x + 1, y) - matrix(x, y)) / mesh.GetStepX();

    const double tdy = (matrix(x, y) - matrix(x, y - 1)) / mesh.GetStepY();
    const double bdy = (matrix(x, y + 1) - matrix(x, y)) / mesh.GetStepY();

    double dx = (ldx - rdx) / mesh.GetStepX();
    double dy = (tdy - bdy) / mesh.GetStepY();

    return dx + dy;
}

void CalculateR(const CMatrix &p, const CMesh &mesh, CMatrix &r)
{
    for( int x = 1; x < r.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < r.GetSizeY() - 1; y++ ) {
            r(x, y) = LaplasOperator(p, mesh, x, y) - RightSide(mesh.X(x), mesh.Y(y));
        }
    }
}

void CalculateG(const CMatrix &r, double alpha, CMatrix &g)
{
    for( int x = 1; x < g.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < g.GetSizeY() - 1; y++ ) {
            g(x, y) = r(x, y) - alpha * g(x, y);
        }
    }
}

double CalculateP(const CMatrix &g, double tau, CMatrix &p, double x_step, double y_step)
{
    double difference = 0;
    for( int x = 1; x < p.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < p.GetSizeY() - 1; y++ ) {
            difference += tau * g(x, y) * tau * g(x, y) * x_step * y_step;
            p(x, y) = p(x, y) - tau * g(x, y);
        }
    }
    return difference;
}

double CalculateAlpha(const CMatrix &r, const CMatrix &g, const CMesh &mesh)
{
    double top = 0., bot = 0.;
    for( int x = 1; x < r.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < r.GetSizeY() - 1; y++ ) {
            double temp = g(x, y) * mesh.GetStepX() * mesh.GetStepY();
            top += LaplasOperator(r, mesh, x, y) * temp;
            bot += LaplasOperator(g, mesh, x, y) * temp;
        }
    }
    return top / bot;
}

double CalculateTau(const CMatrix &r, const CMatrix &g, const CMesh &mesh)
{
    double top = 0., bot = 0.;
    for( int x = 1; x < r.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < r.GetSizeY() - 1; y++ ) {
            double temp = g(x, y) * mesh.GetStepX() * mesh.GetStepY();
            top += r(x, y) * temp;
            bot += LaplasOperator(g, mesh, x, y) * temp;
        }
    }
    return top / bot;
}

double WorkError(const CMatrix &p, const CMesh mesh)
{
    double error = 0;
    for( int x = 1; x < p.GetSizeX() - 1; x++ ) {
        for( int y = 1; y < p.GetSizeY() - 1; y++ ) {
            error = std::max(error, fabs(Boundary(mesh.X(x), mesh.Y(y)) - p(x, y)));
        }
    }
    return error;
}
//

//----CLASSES----
CMatrix::CMatrix(int x_size, int y_size)
{
    m_iSizeX = x_size;
    m_iSizeY = y_size;
    m_Values.resize(m_iSizeX * m_iSizeY);
    std::fill(m_Values.begin(), m_Values.end(), 0.);
}
void zero_elem(double &v)
{
    v = 0;
}
void CMatrix::ClearValues()
{
    std::for_each(m_Values.begin(), m_Values.end(), zero_elem);
}

CMesh::CMesh()
{
    m_iRank = 0;
    m_iProcNumber = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &m_iRank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_iProcNumber);
    m_iGlobalSizeX = 0;
    m_iGlobalSizeY = 0;
    m_dblStepX = 0.;
    m_dblStepY = 0.;
}

void CMesh::ComputeRange(const int point_count, const int block_count, const int index, int &begin, int &end)
{
    int point_per_block = point_count / block_count, rest_points = point_count % block_count;
    begin = index * point_per_block + std::min(index, rest_points);
    end = begin + point_per_block + (index < rest_points ? 1 : 0);
}

long CMesh::Init(const int x_size, const int y_size, const Rectangle all_area)
{
    if( x_size <= 0 || y_size <= 0 )
        return Error;

    m_iGlobalSizeX = x_size;
    m_iGlobalSizeY = y_size;
    m_dblStepX = -(all_area.m_LTPoint.m_dblCrdX - all_area.m_RBPoint.m_dblCrdX) / (double)(x_size - 1);
    m_dblStepY = -(all_area.m_LTPoint.m_dblCrdY - all_area.m_RBPoint.m_dblCrdY) / (double)(y_size - 1);

    if( m_dblStepX <= 0 || m_dblStepY <= 0 )
        return Error;

    int power = 0, temp = 1;
    while( temp < m_iProcNumber ) {
        temp *= 2;
        power++;
    }
    if( temp != m_iProcNumber )
        return Error;

    m_iBlockCountY = pow(2, power / 2);
    if( power%2 ) {
        m_iBlockCountX = pow(2, (power / 2) + 1);
    } else
        m_iBlockCountX = pow(2, power / 2);

    m_iBlockX = m_iRank%m_iBlockCountX;
    m_iBlockY = m_iRank/m_iBlockCountX;

    ComputeRange(x_size, m_iBlockCountX, m_iBlockX, m_iBegX, m_iEndX);
    ComputeRange(y_size, m_iBlockCountY, m_iBlockY, m_iBegY, m_iEndY);
    InitNeibs();

    //std::cout << m_vX.size() << " " << m_vY.size() << std::endl;
    //std::cout << m_iPointPerBlockX * m_iPointPerBlockY << std::endl;

    if( IsNeib(Left) ) {
        m_iBegX--;
    }
    if( IsNeib(Right) ) {
        m_iEndX++;
    }
    if( IsNeib(Top) ) {
        m_iBegY--;
    }
    if( IsNeib(Bot) ) {
        m_iEndY++;
    }

        std::cout << m_iRank << " Block per Y = " << m_iBlockCountY << " Block per X = " << m_iBlockCountX << std::endl;
        std::cout << m_iRank << " P per Y = " << m_iEndY - m_iBegY << " P per X = " << m_iEndX - m_iBegX << std::endl;
        std::cout << m_iRank << " iBegX = " << m_iBegX << " iEndX = " << m_iEndX << std::endl;
        std::cout << m_iRank << " iBegY = " << m_iBegY << " iEndY = " << m_iEndY << std::endl;

    m_vX.clear();
    m_vY.clear();
    m_vX.resize(m_iEndX - m_iBegX);
    m_vY.resize(m_iEndY - m_iBegY);

    //std::cout << m_vX.size() << " " << m_vY.size() << std::endl;
    //std::cout << m_iPointPerBlockX * m_iPointPerBlockY << std::endl;

    for( int i = m_iBegX; i < m_iEndX; i++ ) {
        m_vX[i - m_iBegX] = all_area.m_LTPoint.m_dblCrdX + m_dblStepX * i;
    }

    for( int i = m_iBegY; i < m_iEndY; i++ ) {
        m_vY[i - m_iBegY] = all_area.m_LTPoint.m_dblCrdY + m_dblStepY * i;
    }
    //std::cout << "I "<< m_iRank << " have  Y.l = " << m_vY[0] << " Y.r = " << m_vY[m_vY.size() - 1] << " have  X.t = " << m_vX[0] << " X.b = " << m_vX[m_vX.size() - 1] << std::endl;

    std::vector<MPI_Request> SendRequests, RecieveRequests;
    RecieveRequests.resize(4);
    SendRequests.resize(4);
    std::vector<double> send_buf[4], recieve_buf[4];
    m_vNeibsSizes.resize(4);
    for( int i = 0; i < 4; i++ ) {
        send_buf[i].resize(1);
        recieve_buf[i].resize(1);
        m_vNeibsSizes[i] = 0;
    }

    //std::cout << "I "<< m_iRank << " have  Y = " << m_vY.size() - 2 << " X = " << m_vX.size() - 2 << std::endl;
    for( int i = 0; i < m_vNeibs.size(); i++ ) {
        if( m_vNeibs[i] != -1 ) {
            if( i == (int)Left || i == (int)Right ) {
                send_buf[i][0] = m_vY.size() - 2;
            } else {
                send_buf[i][0] = m_vX.size() - 2;
            }
            //std::cout << "I "<< m_iRank << " send " << send_buf[i][0] << " size\n";
            MPI_Isend(send_buf[i].data(), send_buf[i].size(), MPI_DOUBLE, m_vNeibs[i], 0, MPI_COMM_WORLD, &SendRequests[i]);
            MPI_Irecv(recieve_buf[i].data(), recieve_buf[i].size(), MPI_DOUBLE, m_vNeibs[i], 0, MPI_COMM_WORLD, &RecieveRequests[i]);
        }
    }

    //std::cout << "I "<< m_iRank << " begin w8 for it\n";
    // run wait operations
    for( int i = 0; i < m_vNeibs.size(); i++ ) {
        if( m_vNeibs[i] != -1 ) {
            MPI_Wait(&RecieveRequests[i], MPI_STATUS_IGNORE);
            //std::cout << "I "<< m_iRank << " rec " << recieve_buf[i][0] << " data from" << i << "\n";
            m_vNeibsSizes[i] = recieve_buf[i][0];

            MPI_Wait(&SendRequests[i], MPI_STATUS_IGNORE);
        }
    }

    return OK;
}

void CMesh::InitNeibs()
{
    m_vNeibs.clear();
    m_vNeibs.reserve(4);
    // left nb
    if( m_iBlockX > 0 ) {
        m_vNeibs.push_back((m_iBlockX - 1) + m_iBlockY * m_iBlockCountX);
    } else
        m_vNeibs.push_back(-1);

    // top nb
    if( m_iBlockY > 0 ) {
        m_vNeibs.push_back(m_iBlockX + (m_iBlockY - 1) * m_iBlockCountX);
    } else
        m_vNeibs.push_back(-1);

    // right nb
    if( m_iBlockX < m_iBlockCountX - 1) {
        m_vNeibs.push_back((m_iBlockX + 1) + m_iBlockY * m_iBlockCountX);
    } else
        m_vNeibs.push_back(-1);

    // bot nb
    if( m_iBlockY < m_iBlockCountY - 1) {
        m_vNeibs.push_back(m_iBlockX + (m_iBlockY + 1) * m_iBlockCountX);
    } else
        m_vNeibs.push_back(-1);
}

bool CMesh::IsNeib(Direction dir) const
{
    if( m_vNeibs.size() != 4 )
        return false;

    return ( m_vNeibs[(int)dir] != -1 ? true : false );
}

bool CMesh::Exchange(CMatrix &matrix) const
{
    std::vector<MPI_Request> SendRequests, RecieveRequests;
    RecieveRequests.resize(4);
    SendRequests.resize(4);
    std::vector<double> send_buf[4], recieve_buf[4];
    int size = 0;

    // run asyn operations
    for( int i = 0; i < m_vNeibs.size(); i++ ) {
        if( m_vNeibs[i] != -1 ) {
            recieve_buf[i].resize(m_vNeibsSizes[i]);
            if( i == (int)Left ) {
                size = m_vY.size() - 2;
                send_buf[i].resize(size);
                for(int j = 1; j < m_vY.size() - 1; j++ ) {
                    send_buf[i][j-1] = matrix(1, j);
                }
            } else if( i == (int)Top ) {
                size = m_vX.size() - 2;
                send_buf[i].resize(size);
                for(int j = 1; j < m_vX.size() - 1; j++ ) {
                    send_buf[i][j-1] = matrix(j, 1);
                }
            } else if( i == (int)Right ) {
                size = m_vY.size() - 2;
                send_buf[i].resize(size);
                for(int j = 1; j < m_vY.size() - 1; j++ ) {
                    send_buf[i][j-1] = matrix(m_vX.size() - 2, j);
                }
            } else if( i == (int)Bot ) {
                size = m_vX.size() - 2;
                send_buf[i].resize(size);
                for(int j = 1; j < m_vX.size() - 1; j++ ) {
                    send_buf[i][j-1] = matrix(j, m_vY.size() - 2);
                }
            }
            if( m_iRank == -1 ) {
                std::cout << "I "<< m_iRank << " send " << send_buf[i].size() << " data to " << i;
                for( int k = 0; k < send_buf[i].size(); k++ ) {
                    std::cout << send_buf[i][k] << " ";
                }
                std::cout << std::endl;
            }
            //std::cout << "I "<< m_iRank << " send " << send_buf[i].size() << " data\n";
            MPI_Isend(send_buf[i].data(), send_buf[i].size(), MPI_DOUBLE, m_vNeibs[i], 0, MPI_COMM_WORLD, &SendRequests[i]);
            MPI_Irecv(recieve_buf[i].data(), recieve_buf[i].size(), MPI_DOUBLE, m_vNeibs[i], 0, MPI_COMM_WORLD, &RecieveRequests[i]);
        }
    }

    //std::cout << "I "<< m_iRank << " begin w8 for it\n";
    // run wait operations
    for( int i = 0; i < m_vNeibs.size(); i++ ) {
        if( m_vNeibs[i] != -1 ) {
            MPI_Wait(&RecieveRequests[i], MPI_STATUS_IGNORE);
            //std::cout << "I "<< m_iRank << " rec " << recieve_buf[i].size() << " data\n";
            for( int j = 0; j < recieve_buf[i].size(); j++ ) {
                if( i == (int)Left ) {
                    matrix(0, j + 1) = recieve_buf[i][j];
                } else if( i == (int)Top ) {
                    matrix(j + 1, 0) = recieve_buf[i][j];
                } else if( i == (int)Right ) {
                    matrix(m_vX.size() - 1, j + 1) = recieve_buf[i][j];
                } else if( i == (int)Bot ) {
                    matrix(j + 1, m_vY.size() - 1) = recieve_buf[i][j];
                }
            }

            MPI_Wait(&SendRequests[i], MPI_STATUS_IGNORE);
        }
    }

    return true;
}
//

eReturnCode MainFunction(const int x_size, const int y_size, const Rectangle all_area)
{
    CMesh local_mesh;
    if( local_mesh.Init(x_size, y_size, all_area) != OK )
        return Error;

    //std::cout << local_mesh.GetRank() << " Init ok\n";
    CMatrix p(local_mesh.GetSizeX(), local_mesh.GetSizeY()),
            r(local_mesh.GetSizeX(), local_mesh.GetSizeY()),
            g(local_mesh.GetSizeX(), local_mesh.GetSizeY());

    // init boundary condition
    if( !local_mesh.IsNeib(Left) )
        for( int i = 0; i < p.GetSizeY(); i++ )
            p(0, i) = Boundary(local_mesh.X(0), local_mesh.Y(i));

    if( !local_mesh.IsNeib(Top) )
        for( int i = 0; i < p.GetSizeX(); i++ )
            p(i, 0) = Boundary(local_mesh.X(i), local_mesh.Y(0));

    if( !local_mesh.IsNeib(Right) )
        for( int i = 0; i < p.GetSizeY(); i++ )
            p(p.GetSizeX() - 1, i) = Boundary(local_mesh.X(p.GetSizeX() - 1), local_mesh.Y(i));

    if( !local_mesh.IsNeib(Bot) )
        for( int i = 0; i < p.GetSizeX(); i++ )
            p(i, p.GetSizeY() - 1) = Boundary(local_mesh.X(i), local_mesh.Y(p.GetSizeY() - 1));

   // std::cout << "Start iteration\n";
    double difference = 100, tau = 0., alpha = 0.;
    int iteration = 1;
    while( difference > EPS && difference < 1e10 ) {
        if( local_mesh.GetRank() == 0 ) {
            std::cout << "Iteration #" << iteration++ << std::endl;
            std::cout << "Diff after last iteration = " << difference << std::endl;
        }
        CalculateR(p, local_mesh, r);
        local_mesh.Exchange(r);

        tau = CalculateTau(r, r, local_mesh);
        MPI_Allreduce(MPI_IN_PLACE, &tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        difference = CalculateP(r, tau, p, local_mesh.GetStepX(), local_mesh.GetStepY());
        MPI_Allreduce(MPI_IN_PLACE, &difference, 1 , MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        g = r;
        r.ClearValues();
    }
    while( difference > EPS && difference < 1e10 ) {
        if( local_mesh.GetRank() == 0 ) {
            std::cout << "Iteration #" << iteration++ << std::endl;
            std::cout << "Diff after last iteration = " << difference << std::endl;
        }
        local_mesh.Exchange(p);

        CalculateR(p, local_mesh, r);
        local_mesh.Exchange(r);

        alpha = CalculateAlpha(r, g, local_mesh);
        MPI_Allreduce(MPI_IN_PLACE, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        CalculateG(r, alpha, g);
        local_mesh.Exchange(g);

        tau = CalculateTau(r, g, local_mesh);
        MPI_Allreduce(MPI_IN_PLACE, &tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        difference = CalculateP(g, tau, p, local_mesh.GetStepX(), local_mesh.GetStepY());
        MPI_Allreduce(MPI_IN_PLACE, &difference, 1 , MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    if( local_mesh.GetRank() == 0 ) {
        std::cout << "Max Error = " << WorkError(p, local_mesh) << std::endl;
    }
    if( difference > 1e9 ) {
        if( local_mesh.GetRank() == 0 ) {
            std::cout << "FAIL\n";
        }
        return OK;
    }
    char str[127];
    sprintf(str, "Puasson_CGM_%dx%d.dat", x_size, y_size);
    FILE *fp = fopen(str, "w");

    for (int j = 0; j < x_size; j++)
    {
        for (int i = 0; i < y_size; i++)
            fprintf(fp, "%f %f %f\n", local_mesh.X(j) , local_mesh.Y(i), p(j, i));
    }
    fclose(fp);
    return OK;
}
