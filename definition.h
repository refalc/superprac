#ifndef DEFINITION_H
#define DEFINITION_H
#include <vector>
#include <mpi.h>

const double EPS = 1e-7;
enum eReturnCode
{
    OK = 0,
    Error
};

enum Direction
{
    Left = 0,
    Top,
    Right,
    Bot
};

struct Point2D
{
    double m_dblCrdX;
    double m_dblCrdY;
};

struct Rectangle
{
    Point2D m_LTPoint;         //LeftTopPoint
    Point2D m_RBPoint;         //RightBottomPoint
};

struct IndexInterval
{
    int m_iBeg;
    int m_iEnd;
};

class CMatrix
{
public:
    CMatrix(int x_size, int y_size);
    void ClearValues();
    double& operator()(int x, int y) {return m_Values[x + y * m_iSizeX];}
    double operator()(int x, int y) const {return m_Values[x + y * m_iSizeX];}
    inline int GetSizeX() const {return m_iSizeX;}
    inline int GetSizeY() const {return m_iSizeY;}
private:
    int m_iSizeX;
    int m_iSizeY;
    std::vector<double> m_Values;
};

class CMesh
{
public:
    CMesh();
    long Init(const int x_size, const int y_size, const Rectangle all_area);
    bool IsNeib(Direction dir) const;
    inline double X(int i) const {return m_vX[i];}
    inline double Y(int i) const {return m_vY[i];}
    inline int GetSizeX() const {return m_vX.size();}
    inline int GetSizeY() const {return m_vY.size();}
    inline double GetStepX() const {return m_dblStepX;}
    inline double GetStepY() const {return m_dblStepY;}
    inline int GetRank() const {return m_iRank;}
    bool Exchange(CMatrix &matrix) const;
    void ComputeRange(const int point_count, const int block_count, const int index, int &begin, int &end);
private:
    void InitNeibs();

private:
    int m_iRank;
    int m_iProcNumber;
    int m_iGlobalSizeX;
    int m_iGlobalSizeY;
    int m_iBlockCountX;
    int m_iBlockCountY;
    double m_dblStepX;
    double m_dblStepY;
    IndexInterval m_IntervalX;
    IndexInterval m_IntervalY;
    int m_iBlockX;
    int m_iBlockY;
    int m_iBegX;
    int m_iBegY;
    int m_iEndX;
    int m_iEndY;
    std::vector<int> m_vNeibs;
    std::vector<int> m_vNeibsSizes;
    std::vector<double> m_vX;
    std::vector<double> m_vY;
};

// math functions
double Boundary(double x, double y);
double RightSide(double x, double y);

double LaplasOperator(const CMatrix &matrix, const CMesh &mesh, int x, int y);
void CalculateR(const CMatrix &p, const CMesh &mesh, CMatrix &r);
void CalculateG(const CMatrix &r, double alpha, CMatrix &g);
double CalculateP(const CMatrix &g, double tau, CMatrix &p, double x_step, double y_step);
double CalculateAlpha(const CMatrix &r, const CMatrix &g, const CMesh &mesh);
double CalculateTau(const CMatrix &r, const CMatrix &g, const CMesh &mesh);
double WorkError(const CMatrix &p, const CMesh mesh);

eReturnCode MainFunction(const int x_size, const int y_size, const Rectangle all_area);

#endif // DEFINITION_H
