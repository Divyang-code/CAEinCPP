#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <tuple>
#include <functional>

using namespace std;
const double Pi = 3.141592653589793238;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// This function will set 0.0 in each element of the given vector
void InitVec(int VecSize, double *Vec1, double *Vec2, double *Vec3)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  This function takes 3 vectors (here Vec1, Vec2, and Vec3) and 
    //        their size (here VecSize) as an input.
    // ---->  You may wonder why 3 vectors? Why not one?
    //        Because most of the time in this code, three components (x, y and z)
    //        are need to be initialized. So I wanted to make things abstract.
    //
    // (O-O)  Output:-
    // ---->  3 vectors having 0.0 in each element.

    for(int i = 0; i < VecSize; i++)
    { 
        Vec1[i] = 0.0;
        Vec2[i] = 0.0;
        Vec3[i] = 0.0;
    }
}

// This function will add two vectors of same size
void AddVecs(int Size, double *A, double *B)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  Two vectors (here A and B) of same size (here Size).
    //
    // (O-O)  Output:-
    // ---->  The addition of A and B is stored in the vector A.
    //        That means A is an output of this function.
    //
    // (O=O)  Importance:-
    // ---->  Dummy magentic field in each section i.e. filament, slab, arc, etc.,
    //        we create dummy variable which stores magnetic field value in it
    //        but outside that perticular block of code that dummy variable is unknown
    //        so we add that dummy variable with magnetic field variable which is defined in main().

    for(int i = 0; i < Size; i++)
    {
        A[i] = A[i] + B[i];
    }
}

// This function converts cartesian coordinates into the polar coordinates
tuple<double, double> Cart2Pol(double x, double y)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  Two variables x and y in cartesian coordinates.
    //
    // (O-O)  Output:-
    // ---->  Polar coordinates r and t as tuple.

    double r = pow( pow(x, 2) + pow(y, 2), 0.5);
    double t = atan2(y, x);

    return make_tuple(r, t);
}

// This function converts cartesian mangetic field into the polar field
tuple<double, double> PolarMag(double t, double Bx, double By)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  Cartesian magnetic field's Bx and By components
    //        at a given point which has angle (t) with respect to x-axis.
    //
    // (O-O)  Output:-
    // ---->  Magnetic field in polar fashion.

    double Br = Bx* cos(t) + By* sin(t);
    double Bt =-Bx* sin(t) + By* cos(t);
    
    return make_tuple(Br, Bt);
}

// This functional intigrates the given lambda function
double Trapz(const function<double(double)> &f, int n, double a, double b)
{
    // Functional Details:=
    //
    // (O-O)  Input:-
    // ---->  Lambda function f, number of steps n, starting point a,
    //        and ending point b.
    //
    // (O-O)  Output:-
    // ---->  Intigration of the given lambda function.

    double h = (b - a)/ n;          // Step size

    double SumRest = 0;             // Sum excluding first and last terms

    for(int i = 1; i < n; i++)      // i is running from 1 to (n-1)
    {
        SumRest += f(a + i*h);
    }

    // Trapezoidal rule
    double I = 0.5* h* (f(a) + f(a + n* h) + 2*SumRest);

    return I;
}

// Double integration with trapezoidal rule
double Trapz2D(const function<double(double, double)> &f, int nx, double ax, double bx, int ny, double ay, double by)
{
    // Functional Details:=
    //
    // (O-O)  Input:-
    // ---->  function f with two inputs, number of steps in both x and y direction,
    //        starting and ending values.
    //
    // (O-O)  Output:-
    // ---->  Numerical integration of the given function.

    double hx = (bx - ax)/ nx;                      // step size in x direction
    double hy = (by - ay)/ ny;                      // step size in y direction

    double SumX1 = 0, SumX2 = 0;
    for (int ix = 1; ix < nx; ix++)
    {
        SumX1 += f(ax + ix* hx, ay);
        SumX2 += f(ax + ix* hx, ay + ny* hy);
    }

    double SumY1 = 0, SumY2 = 0;
    for (int iy = 1; iy < ny; iy++)
    {
        SumY1 += f(ax, ay + iy* hy);
        SumY2 += f(ax + nx* hx, ay + iy* hy);
    }

    // x and y summation
    double SumXY = 0;
    for (int ix = 1; ix < nx; ix++)
    {
        for (int iy = 1; iy < ny; iy++)
        {
            SumXY += f(ax + ix* hx, ay + iy* hy);
        }
    }

    // Trapezoidal rule
    double I = (hx* hy/4)* (   f(ax, ay) 
                             + f(ax, ay + ny*hy) 
                             + f(ax + nx*hx, ay) 
                             + f(ax + nx*hx, ay + ny*hy) 
                             + 2*(SumX1 + SumX2 + SumY1 + SumY2) 
                             + 4* SumXY );

    return I;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// This function will load the points where we want to calculate magnetic field
void LoadFieldPoints(int NoOfFieldPoints, 
                     double *Xvector, double *Yvector, double *Zvector)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  This function takes 3 vectors (here Xvector, Yvector, and Zvector),
    //        and their size (here NoOfFieldPoints) as an input.
    //
    // (O-O)  Output:-
    // ---->  At the end of the process, the function will read points.txt text file in which
    //        the coordinates of field points (points where we want to calculate the magnetic field)
    //        is stored.
    // ---->  The x-components of field points are stored in Xvector. Similarly, y and z componets
    //        are stored in Yvector and Zvector, respectively.

    double points[3][NoOfFieldPoints];                  // dummy matrix just to read data

    fstream file("Bridge/points.txt", ios_base::in);    // open the points.txt file (in Bridge folder)

    for (int i = 0; i < 3; i++){                        // variable i is for 3 rows (x, y, and z)
        for (int j = 0; j < NoOfFieldPoints; j++){      // variable j is for colomns or say number of points
            file >> points[i][j];                       // now data is stored in the "points" variable
        }
    }
    file.close();                                       // close the file

    for (int i = 0; i < NoOfFieldPoints; i++){
        Xvector[i] = points[0][i];                      // First row is for x component
        Yvector[i] = points[1][i];                      // Second row is for y component
        Zvector[i] = points[2][i];                      // Third row is for z component
    }
}

// This function will write data to the MAG_DATA.csv file
void WriteDataToCSV(int NoOfFieldPoints, 
                    double *PtX, double *PtY, double *PtZ,
                    double *PtR, double *PtTheta,
                    double *Bx,  double *By,  double *Bz)
{
    // Function Details:=
    //
    // (O-O)  Input:-
    // ---->  This function takes field point and magnetic field vectors, and their size.
    //
    // (O-O)  Output:-
    // ---->  At the end of the process this function will generate MAG_DATA.csv file having
    //        field points and magnetic field data in both cartesian and cylindrical fashion.
    //        In addition, total magnitude of magnetic field is also there.
    //
    // (O=O)  Importance:-
    // ---->  This function plays a vital role as it is not just the output of this function, it is
    //        an output of an entire CAE code.
    // ---->  For further analysis, please open the MAG_DATA.csv file.

    double dummyBr, dummyBt, dummyB;                                                    // dummy variables for radial field Br, angular field Bt, and magnitude B
    
    ofstream MyCSVFile("MAG_DATA.csv");                                                 // Write data to MAG_DATA.csv file

    MyCSVFile << "x (m)" << "\t" << "y (m)" << "\t" << "z (m)" << "\t"                  // first three colomns are for x, y, and z components of field points
              << " " << "\t"                                                            // blank colomn
              << "R (m)" << "\t" << "Theta (deg)" << "\t" << "z (m)" << "\t"            // field points in polar (cylindrical) fashion
              << " " << "\t"
              << "Bx (T)" << "\t" << "By (T)" << "\t" << "Bz (T)" << "\t"               // calculated cartesian magnetic field
              << " " << "\t"
              << "Br (T)" << "\t" << "Btheta (T)" << "\t" << "Bz (T)" << "\t"           // magnetic field in polar (cylindrical) fashion
              << " " << "\t"
              << "B (T)" << "\n";                                                       // total magnetic field - magnitude

    for(int i = 0; i < NoOfFieldPoints; i++)
    {
        tie(dummyBr, dummyBt) = PolarMag(PtTheta[i], Bx[i], By[i]);                     // call PolarMag to find the Br and Bt
        dummyB = pow( pow(Bx[i], 2) + pow(By[i], 2) + pow(Bz[i], 2), 0.5);              // magnitude of total magnetic field

        MyCSVFile << PtX[i] << "\t" << PtY[i] << "\t" << PtZ[i] << "\t"
                  << " " << "\t"
                  << PtR[i] << "\t" << PtTheta[i]* 180/ Pi << "\t" << PtZ[i] << "\t"    // polar points are in degree
                  << " " << "\t"
                  << Bx[i] << "\t" << By[i] << "\t" << Bz[i] << "\t"
                  << " " << "\t"
                  << dummyBr << "\t" << dummyBt << "\t" << Bz[i] << "\t"
                  << " " << "\t"
                  << dummyB << "\n";
    }
    MyCSVFile.close();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Load Filament data
void LoadFilamentData(int NoOfFilament,
                      double *LengthVec, double *IVec,
                      double *XoVec,     double *YoVec,    double *ZoVec,
                      double *OpM11Vec,  double *OpM12Vec, double *OpM13Vec,
                      double *OpM21Vec,  double *OpM22Vec, double *OpM23Vec,
                      double *OpM31Vec,  double *OpM32Vec, double *OpM33Vec)
{
    // Dummy variables to read the data from the file
    double DummyFilVar_center[NoOfFilament][3], DummyFilVar_operator[NoOfFilament][3][3];
    
    fstream filef1("Bridge/f_l.txt",        ios_base::in);
    fstream filef2("Bridge/f_I.txt",        ios_base::in);
    fstream filef3("Bridge/f_center.txt",   ios_base::in);
    fstream filef4("Bridge/f_operator.txt", ios_base::in);
    
    // Read data from the files
    for (int i = 0; i < NoOfFilament; i++)
    {
        filef1 >> LengthVec[i];
        filef2 >> IVec[i];

        for (int j = 0; j < 3; j++){
            filef3 >> DummyFilVar_center[i][j];

            for (int k = 0; k < 3; k++){
                filef4 >> DummyFilVar_operator[i][j][k];
            }
        }
    }
    filef1.close();
    filef2.close();
    filef3.close();
    filef4.close();
    
    // Center and Opeartor variables are initialized
    for (int i = 0; i < NoOfFilament; i++)
    {
        XoVec[i]    = DummyFilVar_center[i][0];
        YoVec[i]    = DummyFilVar_center[i][1];
        ZoVec[i]    = DummyFilVar_center[i][2];

        OpM11Vec[i] = DummyFilVar_operator[i][0][0];
        OpM12Vec[i] = DummyFilVar_operator[i][0][1];
        OpM13Vec[i] = DummyFilVar_operator[i][0][2];
        OpM21Vec[i] = DummyFilVar_operator[i][1][0];
        OpM22Vec[i] = DummyFilVar_operator[i][1][1];
        OpM23Vec[i] = DummyFilVar_operator[i][1][2];
        OpM31Vec[i] = DummyFilVar_operator[i][2][0];
        OpM32Vec[i] = DummyFilVar_operator[i][2][1];
        OpM33Vec[i] = DummyFilVar_operator[i][2][2];
    }

}

// Straight Filament: Mangetic field calculation
class CAEFilamentField
{
private:

    int m_NoOfFilament;                               // Number of filaments (m_ means member)
    int m_NoOfFieldPoint;                             // Number of field data points
    double m_Mu;                                      // permiability

    // Magnetic field calculation for one data point due to one filament only. One to one map.
    tuple<double, double, double> FilamentMagnetic(double X, double Y, double Z,                // coordinates of field point
                                                   double Length, double I,                     // length and current in filament
                                                   double Xo, double Yo, double Zo,             // filament's center
                                                   double OpM11, double OpM12, double OpM13,    // operator variable
                                                   double OpM21, double OpM22, double OpM23,
                                                   double OpM31, double OpM32, double OpM33)
    {
        // Inverse transformation
        double FieldPointX = OpM11* (X-Xo) + OpM21* (Y-Yo) + OpM31* (Z-Zo);
        double FieldPointY = OpM12* (X-Xo) + OpM22* (Y-Yo) + OpM32* (Z-Zo);
        double FieldPointZ = OpM13* (X-Xo) + OpM23* (Y-Yo) + OpM33* (Z-Zo);

        // (xo,yo,z1) and (xo,yo,z2) here xo = 0 and yo = 0
        double zo1 =-Length/2;
        double zo2 = Length/2;
        double xo  = 0.0;
        double yo  = 0.0;

        double u   = FieldPointX - xo;
        double v   = FieldPointY - yo;
        double w1  = FieldPointZ - zo1;
        double w2  = FieldPointZ - zo2;
        double auv = sqrt( pow(u,2) + pow(v,2) );

        double Hx = (v/pow(auv,2))* (-w1/ sqrt( pow(auv,2) + pow(w1,2) ) + w2/ sqrt( pow(auv,2) + pow(w2,2) ) );
        double Hy =-(u/pow(auv,2))* (-w1/ sqrt( pow(auv,2) + pow(w1,2) ) + w2/ sqrt( pow(auv,2) + pow(w2,2) ) );
        double Hz = 0.0;

        double Bix = (m_Mu* I/ (4* Pi))* Hx;
        double Biy = (m_Mu* I/ (4* Pi))* Hy;
        double Biz = (m_Mu* I/ (4* Pi))* Hz;

        // Reverse transformation
        double Bx = OpM11* Bix + OpM12* Biy + OpM13* Biz;
        double By = OpM21* Bix + OpM22* Biy + OpM23* Biz;
        double Bz = OpM31* Bix + OpM32* Biy + OpM33* Biz;
        
        return make_tuple(Bx, By, Bz);
    }

    // Multiple filaments and one data point. Many to one map.
    tuple<double, double, double> ManyFilament(double X, double Y, double Z,
                                               double *LengthVec, double *IVec, 
                                               double *XoVec, double *YoVec, double *ZoVec,
                                               double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                                               double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                                               double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        double Bx = 0.0, By = 0.0, Bz = 0.0;    // Field due to multiple filaments
        double bx, by, bz;                      // Dummy variables

        for(int i = 0; i < m_NoOfFilament; i++)
        {
            tie(bx, by, bz) = FilamentMagnetic(X, Y, Z,                               // remember we need only one data point
                                               LengthVec[i], IVec[i], 
                                               XoVec[i], YoVec[i], ZoVec[i],
                                               OpM11Vec[i], OpM12Vec[i], OpM13Vec[i],
                                               OpM21Vec[i], OpM22Vec[i], OpM23Vec[i],
                                               OpM31Vec[i], OpM32Vec[i], OpM33Vec[i]);

            Bx += bx;
            By += by;
            Bz += bz;
        }

        return make_tuple(Bx, By, Bz);
    }

public:
    
    // Constructor - Initialization.
    CAEFilamentField(int noOfFilament, int noOfFieldPoint, double mu)
    {
        m_NoOfFilament   = noOfFilament;
        m_NoOfFieldPoint = noOfFieldPoint;
        m_Mu             = mu;
    }
    
    // Magnetic field at each data point
    void FilamentMagneticField(double *BxAll, double *ByAll, double *BzAll,
                               double *X, double *Y, double *Z,
                               double *LengthVec, double *IVec, 
                               double *XoVec, double *YoVec, double *ZoVec,
                               double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                               double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                               double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        cout << "Calculating magnetic field due to filaments." << endl;

        for(int i = 0; i < m_NoOfFieldPoint; i++)
        {
            tie(BxAll[i], ByAll[i], BzAll[i]) = ManyFilament(X[i], Y[i], Z[i],
                                                             LengthVec, IVec, 
                                                             XoVec, YoVec, ZoVec,
                                                             OpM11Vec, OpM12Vec, OpM13Vec,
                                                             OpM21Vec, OpM22Vec, OpM23Vec,
                                                             OpM31Vec, OpM32Vec, OpM33Vec);

            cout << i + 1 << "   ";
            if( (i+1)%10 == 0 ){cout << "\n";}
        }
        cout << "\n";
    }
 
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Load Slab data
void LoadSlabData(int NoOfSlab,     double *JVec,
                  double *XVerOrg0, double *YVerOrg0, double *ZVerOrg0,
                  double *XVerOrg6, double *YVerOrg6, double *ZVerOrg6,
                  double *XoVec,    double *YoVec,    double *ZoVec,
                  double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                  double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                  double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
{
    // Dummy variables to read the data from the file
    double DummyVerOrg0[NoOfSlab][3], DummyVerOrg6[NoOfSlab][3], DummySlbVar_center[NoOfSlab][3], DummySlbVar_operator[NoOfSlab][3][3];

    fstream filef1("Bridge/s_J.txt",         ios_base::in);
    fstream filef2("Bridge/s_center.txt",    ios_base::in);
    fstream filef3("Bridge/s_ver_org_0.txt", ios_base::in);
    fstream filef4("Bridge/s_ver_org_6.txt", ios_base::in);
    fstream filef5("Bridge/s_operator.txt",  ios_base::in);

    for (int i = 0; i < NoOfSlab; i++)
    {
        filef1 >> JVec[i];

        for (int j = 0; j < 3; j++){
            filef2 >> DummySlbVar_center[i][j];
            filef3 >> DummyVerOrg0[i][j];
            filef4 >> DummyVerOrg6[i][j];

            for (int k = 0; k < 3; k++){
                filef5 >> DummySlbVar_operator[i][j][k];
            }
        }
    }
    filef1.close();
    filef2.close();
    filef3.close();
    filef4.close();
    filef5.close();

    // Center, Ver_Org, and Opeartor variables are initialized
    for (int i = 0; i < NoOfSlab; i++)
    {    
        XoVec[i]    = DummySlbVar_center[i][0];
        YoVec[i]    = DummySlbVar_center[i][1];
        ZoVec[i]    = DummySlbVar_center[i][2];

        XVerOrg0[i] = DummyVerOrg0[i][0];
        YVerOrg0[i] = DummyVerOrg0[i][1];
        ZVerOrg0[i] = DummyVerOrg0[i][2];

        XVerOrg6[i] = DummyVerOrg6[i][0];
        YVerOrg6[i] = DummyVerOrg6[i][1];
        ZVerOrg6[i] = DummyVerOrg6[i][2];

        OpM11Vec[i] = DummySlbVar_operator[i][0][0];
        OpM12Vec[i] = DummySlbVar_operator[i][0][1];
        OpM13Vec[i] = DummySlbVar_operator[i][0][2];
        OpM21Vec[i] = DummySlbVar_operator[i][1][0];
        OpM22Vec[i] = DummySlbVar_operator[i][1][1];
        OpM23Vec[i] = DummySlbVar_operator[i][1][2];
        OpM31Vec[i] = DummySlbVar_operator[i][2][0];
        OpM32Vec[i] = DummySlbVar_operator[i][2][1];
        OpM33Vec[i] = DummySlbVar_operator[i][2][2];
    }

}

// Straight Slab: Magnetic field calculation
class CAESlabField
{
private:
    int m_NoOfSlab;                                   // Number of slabs
    int m_NoOfFieldPoint;                             // Number of field data points
    double m_Mu;                                      // permiability

    double tau(double ui, double vj, double wk)
    {
        double Val = -(  wk* asinh( ui/ sqrt( pow(wk, 2) + pow(vj, 2) ) )
                       + ui* asinh( wk/ sqrt( pow(ui, 2) + pow(vj, 2) ) )
                       - vj* atan( (ui/ vj)* (wk/ sqrt( pow(ui, 2) + pow(vj, 2) + pow(wk, 2)) ) ) 
                      );
    
        return Val;
    }

    // Magnetic field calculation for one data point due to one slab only. One to one map.
    tuple<double, double, double> SlabMagnetic(double X, double Y, double Z,
                                               double XVR0, double YVR0, double ZVR0,
                                               double XVR6, double YVR6, double ZVR6,
                                               double J,
                                               double Xo, double Yo, double Zo,
                                               double OpM11, double OpM12, double OpM13,
                                               double OpM21, double OpM22, double OpM23,
                                               double OpM31, double OpM32, double OpM33)
    {
        // Inverse transformation
        double FieldPointX = OpM11* (X-Xo) + OpM21* (Y-Yo) + OpM31* (Z-Zo);
        double FieldPointY = OpM12* (X-Xo) + OpM22* (Y-Yo) + OpM32* (Z-Zo);
        double FieldPointZ = OpM13* (X-Xo) + OpM23* (Y-Yo) + OpM33* (Z-Zo);

        double u1 = FieldPointX - XVR0;
        double u2 = FieldPointX - XVR6;

        double v1 = FieldPointY - YVR0;
        double v2 = FieldPointY - YVR6;

        double w1 = FieldPointZ - ZVR0;
        double w2 = FieldPointZ - ZVR6;

        double Hx = ( - tau(u1, v1, w1) 
                      + tau(u1, v1, w2) 
                      + tau(u1, v2, w1) 
                      - tau(u1, v2, w2)
                      + tau(u2, v1, w1) 
                      - tau(u2, v1, w2) 
                      - tau(u2, v2, w1) 
                      + tau(u2, v2, w2) );
        
        double Hy = - ( - tau(v1, u1, w1)
                        + tau(v1, u1, w2)
                        + tau(v1, u2, w1)
                        - tau(v1, u2, w2)
                        + tau(v2, u1, w1)
                        - tau(v2, u1, w2)
                        - tau(v2, u2, w1)
                        + tau(v2, u2, w2) );
        
        double Hz = 0;


        double Bix = (m_Mu* J/ (4* Pi))* Hx;
        double Biy = (m_Mu* J/ (4* Pi))* Hy;
        double Biz = (m_Mu* J/ (4* Pi))* Hz;

        // Reverse transformation
        double Bx = OpM11* Bix + OpM12* Biy + OpM13* Biz;
        double By = OpM21* Bix + OpM22* Biy + OpM23* Biz;
        double Bz = OpM31* Bix + OpM32* Biy + OpM33* Biz;

        return make_tuple(Bx, By, Bz);
    }

    // Multiple slabs and one data point. Many to one map.
    tuple<double, double, double> ManySlab(double X, double Y, double Z,
                                           double *XVR0Vec, double *YVR0Vec, double *ZVR0Vec,
                                           double *XVR6Vec, double *YVR6Vec, double *ZVR6Vec,
                                           double *JVec, 
                                           double *XoVec, double *YoVec, double *ZoVec,
                                           double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                                           double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                                           double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        double Bx = 0.0, By = 0.0, Bz = 0.0;    // Field due to multiple slabs
        double bx, by, bz;                      // Dummy variables

        for(int i = 0; i < m_NoOfSlab; i++)
        {
            tie(bx, by, bz) = SlabMagnetic(X, Y, Z,                               // remember we need only one data point
                                           XVR0Vec[i], YVR0Vec[i], ZVR0Vec[i],
                                           XVR6Vec[i], YVR6Vec[i], ZVR6Vec[i],
                                           JVec[i], 
                                           XoVec[i], YoVec[i], ZoVec[i],
                                           OpM11Vec[i], OpM12Vec[i], OpM13Vec[i],
                                           OpM21Vec[i], OpM22Vec[i], OpM23Vec[i],
                                           OpM31Vec[i], OpM32Vec[i], OpM33Vec[i]);

            Bx += bx;
            By += by;
            Bz += bz;
        }

        return make_tuple(Bx, By, Bz);
    }

public:
    
    // Constructor - Initialization.
    CAESlabField(int noOfSlab, int noOfFieldPoint, double mu)
    {
        m_NoOfSlab       = noOfSlab;
        m_NoOfFieldPoint = noOfFieldPoint;
        m_Mu             = mu;
    }
    
    // Magnetic field at each data point
    void SlabMagneticField(double *BxAll, double *ByAll, double *BzAll,
                               double *X, double *Y, double *Z,
                               double *XVR0Vec, double *YVR0Vec, double *ZVR0Vec,
                               double *XVR6Vec, double *YVR6Vec, double *ZVR6Vec,
                               double *JVec,
                               double *XoVec, double *YoVec, double *ZoVec,
                               double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                               double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                               double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        cout << "Calculating magnetic field due to Slabs." << endl;

        for(int i = 0; i < m_NoOfFieldPoint; i++)
        {
            tie(BxAll[i], ByAll[i], BzAll[i]) = ManySlab(X[i], Y[i], Z[i],
                                                         XVR0Vec, YVR0Vec, ZVR0Vec,
                                                         XVR6Vec, YVR6Vec, ZVR6Vec,
                                                         JVec,
                                                         XoVec, YoVec, ZoVec,
                                                         OpM11Vec, OpM12Vec, OpM13Vec,
                                                         OpM21Vec, OpM22Vec, OpM23Vec,
                                                         OpM31Vec, OpM32Vec, OpM33Vec);

            cout << i + 1 << "   ";
            if( (i+1)%10 == 0 ){cout << "\n";}
        }
        cout << "\n";
    }
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Load Circular Arc data
void LoadArcData(int NoOfArc, 
                 double *IVec,     double *a0Vec,
                 double *Phy1Vec,  double *Phy2Vec,
                 double *XoVec,    double *YoVec,    double *ZoVec,
                 double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                 double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                 double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
{
    // Dummy variables to read the data from the file
    double DummyArcVar_center[NoOfArc][3], DummyArcVar_operator[NoOfArc][3][3];

    fstream filef1("Bridge/a_I.txt",        ios_base::in);
    fstream filef2("Bridge/a_cnt_cart.txt", ios_base::in);
    fstream filef3("Bridge/a_phy_1.txt",    ios_base::in);
    fstream filef4("Bridge/a_phy_2.txt",    ios_base::in);
    fstream filef5("Bridge/a_operator.txt", ios_base::in);
    fstream filef6("Bridge/a_a0.txt",       ios_base::in);

    for (int i = 0; i < NoOfArc; i++)
    {
        filef1 >> IVec[i];
        filef3 >> Phy1Vec[i];
        filef4 >> Phy2Vec[i];
        filef6 >> a0Vec[i];

        for (int j = 0; j < 3; j++){
            filef2 >> DummyArcVar_center[i][j];

            for (int k = 0; k < 3; k++){
                filef5 >> DummyArcVar_operator[i][j][k];
            }
        }
    }
    filef1.close();
    filef2.close();
    filef3.close();
    filef4.close();
    filef5.close();
    filef6.close();

    // Center, and Opeartor variables are initialized
    for (int i = 0; i < NoOfArc; i++)
    {    
        XoVec[i]    = DummyArcVar_center[i][0];
        YoVec[i]    = DummyArcVar_center[i][1];
        ZoVec[i]    = DummyArcVar_center[i][2];

        OpM11Vec[i] = DummyArcVar_operator[i][0][0];
        OpM12Vec[i] = DummyArcVar_operator[i][0][1];
        OpM13Vec[i] = DummyArcVar_operator[i][0][2];
        OpM21Vec[i] = DummyArcVar_operator[i][1][0];
        OpM22Vec[i] = DummyArcVar_operator[i][1][1];
        OpM23Vec[i] = DummyArcVar_operator[i][1][2];
        OpM31Vec[i] = DummyArcVar_operator[i][2][0];
        OpM32Vec[i] = DummyArcVar_operator[i][2][1];
        OpM33Vec[i] = DummyArcVar_operator[i][2][2];
    }
}

// Circular Arc: Magnetic field calculation
class CAEArcField
{
private:
    int m_NoOfArc;                                    // Number of Arcs
    int m_NoOfFieldPoint;                             // Number of field data points
    double m_Mu;                                      // permiability

    // Magnetic field calculation for one data point due to one arc only. One to one map.
    tuple<double, double, double> ArcMagnetic(double X, double Y, double Z,                // coordinates of field point
                                              double a0, double I,                         // radius and current in arc
                                              double phy1, double phy2,                    // initial and final angle of arc
                                              double Xo, double Yo, double Zo,             // arc's center
                                              double OpM11, double OpM12, double OpM13,    // operator variable
                                              double OpM21, double OpM22, double OpM23,
                                              double OpM31, double OpM32, double OpM33)
    {
        // Inverse transformation
        double FP_X = OpM11* (X-Xo) + OpM21* (Y-Yo) + OpM31* (Z-Zo);
        double FP_Y = OpM12* (X-Xo) + OpM22* (Y-Yo) + OpM32* (Z-Zo);
        double FP_Z = OpM13* (X-Xo) + OpM23* (Y-Yo) + OpM33* (Z-Zo);

        double FP_R, FP_Phy;                                                      // FP means Field Point
        tie(FP_R, FP_Phy) = Cart2Pol(FP_X, FP_Y);                                 // field point in polar fashion

        // Constants
        double b = sqrt( pow(a0, 2) + pow(FP_R, 2) + pow(FP_Z, 2) );
        double Phi1 = phy1 - FP_Phy;
        double Phi2 = phy2 - FP_Phy;

        // Lambda functions for integration I1 and I3
        auto lmbda_In1 = [b, FP_R, a0](double PhiN){ return cos(PhiN)/ pow( pow(b, 2) - 2* FP_R* a0* cos(PhiN), 1.5); };
        auto lmbda_In3 = [b, FP_R, a0](double PhiN){ return 1/ pow( pow(b, 2) - 2* FP_R* a0* cos(PhiN), 1.5); };

        // To find the value of In1 and In2, integrate them numerically
        double In1 = Trapz(lmbda_In1, 1000, Phi1, Phi2);
        double In3 = Trapz(lmbda_In3, 1000, Phi1, Phi2);
        
        // Value of In2
        double In2;

        if (FP_R == 0){
            In2 = (1/ pow( pow(a0, 2) + pow(FP_Z, 2), 1.5) ) * (cos(Phi1) - cos(Phi2));
        }
        else{
            In2 =-(1/(FP_R* a0))* (1/ sqrt( pow(b, 2) - 2* FP_R* a0* cos(Phi2)) - 1/sqrt( pow(b, 2) - 2* FP_R* a0* cos(Phi1)) );
        }

        double Hx = a0* FP_Z* ( cos(FP_Phy)* In1 - sin(FP_Phy)* In2);
        double Hy = a0* FP_Z* ( sin(FP_Phy)* In1 + cos(FP_Phy)* In2);
        double Hz = a0* (a0* In3 - FP_R* In1);

        double Bix = (m_Mu* I/ (4* Pi))* Hx;
        double Biy = (m_Mu* I/ (4* Pi))* Hy;
        double Biz = (m_Mu* I/ (4* Pi))* Hz;

        // Reverse transformation
        double Bx = OpM11* Bix + OpM12* Biy + OpM13* Biz;
        double By = OpM21* Bix + OpM22* Biy + OpM23* Biz;
        double Bz = OpM31* Bix + OpM32* Biy + OpM33* Biz;

        return make_tuple(Bx, By, Bz);
    }

    // Multiple filaments and one data point. Many to one map.
    tuple<double, double, double> ManyArc(double X, double Y, double Z,
                                          double *a0Vec, double *IVec,
                                          double *phy1Vec, double *phy2Vec,
                                          double *XoVec, double *YoVec, double *ZoVec,
                                          double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                                          double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                                          double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        double Bx = 0.0, By = 0.0, Bz = 0.0;    // Field due to multiple filaments
        double bx, by, bz;                      // Dummy variables

        for(int i = 0; i < m_NoOfArc; i++)
        {
            tie(bx, by, bz) = ArcMagnetic(X, Y, Z,                               // remember we need only one data point
                                          a0Vec[i], IVec[i],
                                          phy1Vec[i], phy2Vec[i],
                                          XoVec[i], YoVec[i], ZoVec[i],
                                          OpM11Vec[i], OpM12Vec[i], OpM13Vec[i],
                                          OpM21Vec[i], OpM22Vec[i], OpM23Vec[i],
                                          OpM31Vec[i], OpM32Vec[i], OpM33Vec[i]);

            Bx += bx;
            By += by;
            Bz += bz;
        }

        return make_tuple(Bx, By, Bz);
    }

public:

    CAEArcField(int noOfArc, int noOfFieldPoint, double mu)
    {
        m_NoOfArc        = noOfArc;
        m_NoOfFieldPoint = noOfFieldPoint;
        m_Mu             = mu;
    }
    
    // Magnetic field at each data point
    void ArcMagneticField(double *BxAll, double *ByAll, double *BzAll,
                          double *X, double *Y, double *Z,
                          double *a0Vec, double *IVec,
                          double *phy1Vec, double *phy2Vec,
                          double *XoVec, double *YoVec, double *ZoVec,
                          double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                          double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                          double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        cout << "Calculating magnetic field due to arcs." << endl;

        for(int i = 0; i < m_NoOfFieldPoint; i++)
        {
            tie(BxAll[i], ByAll[i], BzAll[i]) = ManyArc(X[i], Y[i], Z[i],
                                                        a0Vec, IVec,
                                                        phy1Vec, phy2Vec,
                                                        XoVec, YoVec, ZoVec,
                                                        OpM11Vec, OpM12Vec, OpM13Vec,
                                                        OpM21Vec, OpM22Vec, OpM23Vec,
                                                        OpM31Vec, OpM32Vec, OpM33Vec);

            cout << i + 1 << "   ";
            if( (i+1)%10 == 0 ){cout << "\n";}
        }
        cout << "\n";
    }

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Load Rectangular Arc data
void LoadRectArcData(int NoOfRectArc, 
                     double *RInVec,   double *ROutVec,
                     double *TcknsVec, double *JVec,
                     double *Phy1Vec,  double *Phy2Vec,
                     double *XoVec,    double *YoVec,    double *ZoVec,
                     double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                     double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                     double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
{
    // Dummy variables to read the data from the file
    double DummyRAVar_CntCart[NoOfRectArc][3], DummyRAVar_operator[NoOfRectArc][3][3];

    fstream filera1("Bridge/ra_r_inner.txt",  ios_base::in);
    fstream filera2("Bridge/ra_r_outer.txt",  ios_base::in);
    fstream filera3("Bridge/ra_tckns.txt",    ios_base::in);
    fstream filera4("Bridge/ra_phy_1.txt",    ios_base::in);
    fstream filera5("Bridge/ra_phy_2.txt",    ios_base::in);
    fstream filera6("Bridge/ra_J.txt",        ios_base::in);
    fstream filera7("Bridge/ra_cnt_cart.txt", ios_base::in);
    fstream filera8("Bridge/ra_operator.txt", ios_base::in);

    for (int i = 0; i < NoOfRectArc; i++)
    {
        filera1 >> RInVec[i];
        filera2 >> ROutVec[i];
        filera3 >> TcknsVec[i];
        filera4 >> Phy1Vec[i];
        filera5 >> Phy2Vec[i];
        filera6 >> JVec[i];

        for (int j = 0; j < 3; j++)
        {
            filera7 >> DummyRAVar_CntCart[i][j];

            for (int k = 0; k < 3; k++)
            {
                filera8 >> DummyRAVar_operator[i][j][k];
            }
        }
    }
    filera1.close();
    filera2.close();
    filera3.close();
    filera4.close();
    filera5.close();
    filera6.close();
    filera7.close();
    filera8.close();

    // Center and Opeartor variables are initialized
    for (int i = 0; i < NoOfRectArc; i++)
    {
        XoVec[i]    = DummyRAVar_CntCart[i][0];
        YoVec[i]    = DummyRAVar_CntCart[i][1];
        ZoVec[i]    = DummyRAVar_CntCart[i][2];

        OpM11Vec[i] = DummyRAVar_operator[i][0][0];
        OpM12Vec[i] = DummyRAVar_operator[i][0][1];
        OpM13Vec[i] = DummyRAVar_operator[i][0][2];
        OpM21Vec[i] = DummyRAVar_operator[i][1][0];
        OpM22Vec[i] = DummyRAVar_operator[i][1][1];
        OpM23Vec[i] = DummyRAVar_operator[i][1][2];
        OpM31Vec[i] = DummyRAVar_operator[i][2][0];
        OpM32Vec[i] = DummyRAVar_operator[i][2][1];
        OpM33Vec[i] = DummyRAVar_operator[i][2][2];
    }
}

// Circular Arc with rectangular cross-section
class CAERectArcField
{
private:
    int m_NoOfRectArc;                                // Number of Rectangular Arcs
    int m_NoOfFieldPoint;                             // Number of field data points
    double m_Mu;                                      // permiability
    int m_Phi_Steps;                                  // number of steps for angle phi
    double m_Acc_r;                                   // accuracy for r

    // I_Rr1
    double Solve_IRr1(double FP_R, int Steps,
                      double r1,   double r2, 
                      double w1,   double w2,
                      double Phi1, double Phi2)
    {
        // I_Rr1 : Without integration
        auto lmbda_I_Rr1_11 = [FP_R, r1, w1](double PhiN)
        {  
            return  cos(PhiN)* sqrt( pow(FP_R, 2) + pow(r1, 2) + pow(w1, 2) - 2* FP_R* r1* cos(PhiN) )
                    + FP_R* pow(cos(PhiN), 2)* log( r1 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r1, 2) + pow(w1, 2) - 2* FP_R* r1* cos(PhiN)) ); 
        };

        auto lmbda_I_Rr1_12 = [FP_R, r1, w2](double PhiN)
        {  
            return  cos(PhiN)* sqrt( pow(FP_R, 2) + pow(r1, 2) + pow(w2, 2) - 2* FP_R* r1* cos(PhiN) )
                    + FP_R* pow(cos(PhiN), 2)* log( r1 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r1, 2) + pow(w2, 2) - 2* FP_R* r1* cos(PhiN)) ); 
        };

        auto lmbda_I_Rr1_21 = [FP_R, r2, w1](double PhiN)
        {  
            return  cos(PhiN)* sqrt( pow(FP_R, 2) + pow(r2, 2) + pow(w1, 2) - 2* FP_R* r2* cos(PhiN) )
                    + FP_R* pow(cos(PhiN), 2)* log( r2 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r2, 2) + pow(w1, 2) - 2* FP_R* r2* cos(PhiN)) ); 
        };

        auto lmbda_I_Rr1_22 = [FP_R, r2, w2](double PhiN)
        {  
            return  cos(PhiN)* sqrt( pow(FP_R, 2) + pow(r2, 2) + pow(w2, 2) - 2* FP_R* r2* cos(PhiN) )
                    + FP_R* pow(cos(PhiN), 2)* log( r2 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r2, 2) + pow(w2, 2) - 2* FP_R* r2* cos(PhiN)) ); 
        };

        // I_Rr1 : Integration done
        double I_Rr1_11 = Trapz(lmbda_I_Rr1_11, Steps, Phi1, Phi2);
        double I_Rr1_12 = Trapz(lmbda_I_Rr1_12, Steps, Phi1, Phi2);
        double I_Rr1_21 = Trapz(lmbda_I_Rr1_21, Steps, Phi1, Phi2);
        double I_Rr1_22 = Trapz(lmbda_I_Rr1_22, Steps, Phi1, Phi2);

        double I_Rr1 =   I_Rr1_11
                       - I_Rr1_12
                       - I_Rr1_21
                       + I_Rr1_22;

        return I_Rr1;
    }

    // I_Rr2
    double Solve_IRr2(double FP_R, double ri, double phi_i, double wk)
    {
        double I_Rr2;

        if (FP_R == 0.0)
        {
            I_Rr2 = - cos(phi_i)* sqrt(pow(ri, 2) + pow(wk, 2));
        }
        else
        {
            I_Rr2 = (0.5/FP_R)* ( (ri - FP_R* cos(phi_i))* sqrt(pow(FP_R, 2) + pow(ri, 2) + pow(wk, 2) - 2* FP_R* ri* cos(phi_i))
                    + (pow(wk, 2) + pow(FP_R* sin(phi_i), 2))* log(ri - FP_R* cos(phi_i) + sqrt(pow(FP_R, 2) + pow(ri, 2) + pow(wk, 2) - 2* FP_R* ri* cos(phi_i)) ) );
        }
        
        return I_Rr2;
    }

    // I_Zr1
    double Solve_IZr1(double FP_R, int Steps,
                      double r1,   double r2, 
                      double w1,   double w2,
                      double Phi1, double Phi2)
    {
        // I_Zr1 : Without integration
        auto lmbda_I_Zr1_11 = [FP_R, r1, w1](double PhiN)
        {
            return w1* log(r1 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r1, 2) + pow(w1, 2) - 2* FP_R* r1* cos(PhiN)) );
        };

        auto lmbda_I_Zr1_12 = [FP_R, r1, w2](double PhiN)
        {
            return w2* log(r1 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r1, 2) + pow(w2, 2) - 2* FP_R* r1* cos(PhiN)) );
        };

        auto lmbda_I_Zr1_21 = [FP_R, r2, w1](double PhiN)
        {
            return w1* log(r2 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r2, 2) + pow(w1, 2) - 2* FP_R* r2* cos(PhiN)) );
        };

        auto lmbda_I_Zr1_22 = [FP_R, r2, w2](double PhiN)
        {
            return w2* log(r2 - FP_R* cos(PhiN) + sqrt(pow(FP_R, 2) + pow(r2, 2) + pow(w2, 2) - 2* FP_R* r2* cos(PhiN)) );
        };

        double I_Zr1_11 = Trapz(lmbda_I_Zr1_11, Steps, Phi1, Phi2);
        double I_Zr1_12 = Trapz(lmbda_I_Zr1_12, Steps, Phi1, Phi2);
        double I_Zr1_21 = Trapz(lmbda_I_Zr1_21, Steps, Phi1, Phi2);
        double I_Zr1_22 = Trapz(lmbda_I_Zr1_22, Steps, Phi1, Phi2);

        double I_Zr1 =   I_Zr1_11
                       - I_Zr1_12
                       - I_Zr1_21
                       + I_Zr1_22;

        return I_Zr1;
    }

    // I_Zr2
    double Solve_IZr2(double FP_R, 
                      double AccR, int StepP, 
                      double r1,   double r2, 
                      double w1,   double w2,
                      double Phi1, double Phi2)
    {
        // I_Zr2 : Without integration
        auto lmbda_I_Zr2_1 = [FP_R, w1](double rN, double PhiN)
        {
            return w1* (pow(FP_R, 2) - FP_R* rN* cos(PhiN))/ ( (pow(FP_R, 2) + pow(rN, 2) - 2* FP_R* rN* cos(PhiN))* sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w1, 2) - 2* FP_R* rN* cos(PhiN)) );
        };

        auto lmbda_I_Zr2_2 = [FP_R, w2](double rN, double PhiN)
        {
            return w2* (pow(FP_R, 2) - FP_R* rN* cos(PhiN))/ ( (pow(FP_R, 2) + pow(rN, 2) - 2* FP_R* rN* cos(PhiN))* sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w2, 2) - 2* FP_R* rN* cos(PhiN)) );
        };

        int StepR = 1 + (r2 - r1)/ AccR;

        double I_Zr2_1 = Trapz2D(lmbda_I_Zr2_1, StepR, r1, r2, StepP, Phi1, Phi2);
        double I_Zr2_2 = Trapz2D(lmbda_I_Zr2_2, StepR, r1, r2, StepP, Phi1, Phi2);

        double I_Zr2 = - I_Zr2_1 + I_Zr2_2;

        return I_Zr2;
    }

    // Magnetic field calculation for one data point due to one filament only. One to one map.
    tuple<double, double, double> RectArcMagnetic(double X,     double Y,     double Z,                // coordinates of field point
                                                  double RIn1,  double ROut2,                         // inner and outter and radius
                                                  double Tckns, double J,                              // thickness and current density 
                                                  double Phy1,  double Phy2,                           // initial and final angle of rect. arc 
                                                  double Xo,    double Yo,    double Zo,               // center
                                                  double OpM11, double OpM12, double OpM13,            // operator variable
                                                  double OpM21, double OpM22, double OpM23,
                                                  double OpM31, double OpM32, double OpM33)
    {
        // Inverse transformation
        double FP_X = OpM11* (X-Xo) + OpM21* (Y-Yo) + OpM31* (Z-Zo);        // FP means Field Point
        double FP_Y = OpM12* (X-Xo) + OpM22* (Y-Yo) + OpM32* (Z-Zo);
        double FP_Z = OpM13* (X-Xo) + OpM23* (Y-Yo) + OpM33* (Z-Zo);

        double FP_R, FP_Phy;
        tie(FP_R, FP_Phy) = Cart2Pol(FP_X, FP_Y);                           // field point in polar fashion

        double z1 = -Tckns/2;
        double z2 =  Tckns/2;

        double Phi1 = Phy1 - FP_Phy;
        double Phi2 = Phy2 - FP_Phy;
        double w1   = FP_Z - z1;
        double w2   = FP_Z - z2;

        double I_Rr1 = Solve_IRr1(FP_R, m_Phi_Steps, RIn1, ROut2, w1, w2, Phi1, Phi2);
        double I_Rr2 = - Solve_IRr2(FP_R, RIn1,  Phi1, w1)
                       + Solve_IRr2(FP_R, RIn1,  Phi1, w2)
                       + Solve_IRr2(FP_R, RIn1,  Phi2, w1)
                       - Solve_IRr2(FP_R, RIn1,  Phi2, w2)
                       + Solve_IRr2(FP_R, ROut2, Phi1, w1)
                       - Solve_IRr2(FP_R, ROut2, Phi1, w2)
                       - Solve_IRr2(FP_R, ROut2, Phi2, w1)
                       + Solve_IRr2(FP_R, ROut2, Phi2, w2);
        
        double I_Zr1 = Solve_IZr1(FP_R, m_Phi_Steps, RIn1, ROut2, w1, w2, Phi1, Phi2);
        double I_Zr2 = Solve_IZr2(FP_R, m_Acc_r, m_Phi_Steps, RIn1, ROut2, w1, w2, Phi1, Phi2);

        double Hx = cos(FP_Phy)* I_Rr1 - sin(FP_Phy)* I_Rr2;
        double Hy = cos(FP_Phy)* I_Rr2 + sin(FP_Phy)* I_Rr1;
        double Hz = - I_Zr1 + I_Zr2;

        double Bix = (m_Mu* J/ (4* Pi))* Hx;
        double Biy = (m_Mu* J/ (4* Pi))* Hy;
        double Biz = (m_Mu* J/ (4* Pi))* Hz;

        // Reverse transformation
        double Bx = OpM11* Bix + OpM12* Biy + OpM13* Biz;
        double By = OpM21* Bix + OpM22* Biy + OpM23* Biz;
        double Bz = OpM31* Bix + OpM32* Biy + OpM33* Biz;
        
        return make_tuple(Bx, By, Bz);
    }

    // Multiple filaments and one data point. Many to one map.
    tuple<double, double, double> ManyRectArc(double X,         double Y,         double Z,
                                              double *RInVec,   double *ROutVec,
                                              double *TcknsVec, double *JVec,
                                              double *Phy1Vec,  double *Phy2Vec,
                                              double *XoVec,    double *YoVec,    double *ZoVec,
                                              double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                                              double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                                              double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        double Bx = 0.0, By = 0.0, Bz = 0.0;    // Field due to multiple filaments
        double bx, by, bz;                      // Dummy variables

        for(int i = 0; i < m_NoOfRectArc; i++)
        {
            tie(bx, by, bz) = RectArcMagnetic(X, Y, Z,                                  // remember we need only one data point
                                              RInVec[i],   ROutVec[i],
                                              TcknsVec[i], JVec[i],
                                              Phy1Vec[i],  Phy2Vec[i],
                                              XoVec[i],    YoVec[i],    ZoVec[i],
                                              OpM11Vec[i], OpM12Vec[i], OpM13Vec[i],
                                              OpM21Vec[i], OpM22Vec[i], OpM23Vec[i],
                                              OpM31Vec[i], OpM32Vec[i], OpM33Vec[i]);

            Bx += bx;
            By += by;
            Bz += bz;
        }

        return make_tuple(Bx, By, Bz);
    }

public:

    CAERectArcField(int noOfRectArc, int noOfFieldPoint, double mu, int phi_Steps, double acc_r)
    {
        m_NoOfRectArc    = noOfRectArc;
        m_NoOfFieldPoint = noOfFieldPoint;
        m_Mu             = mu;
        m_Phi_Steps      = phi_Steps;
        m_Acc_r          = acc_r;
    }

    // Magnetic field at each data point
    void RectArcMagneticField(double *BxAll,    double *ByAll,    double *BzAll,
                              double *X,        double *Y,        double *Z,
                              double *RInVec,   double *ROutVec,
                              double *TcknsVec, double *JVec,
                              double *Phy1Vec,  double *Phy2Vec,
                              double *XoVec,    double *YoVec,    double *ZoVec,
                              double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                              double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                              double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        cout << "Calculating magnetic field due to rect. arcs." << endl;

        for(int i = 0; i < m_NoOfFieldPoint; i++)
        {
            tie(BxAll[i], ByAll[i], BzAll[i]) = ManyRectArc(X[i], Y[i], Z[i],
                                                            RInVec,   ROutVec,
                                                            TcknsVec, JVec,
                                                            Phy1Vec,  Phy2Vec,
                                                            XoVec,    YoVec,    ZoVec,
                                                            OpM11Vec, OpM12Vec, OpM13Vec,
                                                            OpM21Vec, OpM22Vec, OpM23Vec,
                                                            OpM31Vec, OpM32Vec, OpM33Vec);

            cout << i + 1 << "   ";
            if( (i+1)%10 == 0 ){cout << "\n";}
        }
        cout << "\n";
    }

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Load Cylinder data
void LoadCyldrData(int NoOfCylndr, 
                   double *RoVec,    double *LVec,     double *JVec,
                   double *XoVec,    double *YoVec,    double *ZoVec,
                   double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                   double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                   double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
{
    // Dummy variables to read the data from the file
    double DummyRAVar_CntCart[NoOfCylndr][3], DummyRAVar_operator[NoOfCylndr][3][3];

    fstream file1("Bridge/cy_ro.txt",       ios_base::in);
    fstream file2("Bridge/cy_l.txt",        ios_base::in);
    fstream file3("Bridge/cy_J.txt",        ios_base::in);
    fstream file4("Bridge/cy_center.txt",   ios_base::in);
    fstream file5("Bridge/cy_operator.txt", ios_base::in);

    for (int i = 0; i < NoOfCylndr; i++)
    {
        file1 >> RoVec[i];
        file2 >> LVec[i];
        file3 >> JVec[i];

        for (int j = 0; j < 3; j++)
        {
            file4 >> DummyRAVar_CntCart[i][j];

            for (int k = 0; k < 3; k++)
            {
                file5 >> DummyRAVar_operator[i][j][k];
            }
        }
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();

    // Center and Opeartor variables are initialized
    for (int i = 0; i < NoOfCylndr; i++)
    {
        XoVec[i]    = DummyRAVar_CntCart[i][0];
        YoVec[i]    = DummyRAVar_CntCart[i][1];
        ZoVec[i]    = DummyRAVar_CntCart[i][2];

        OpM11Vec[i] = DummyRAVar_operator[i][0][0];
        OpM12Vec[i] = DummyRAVar_operator[i][0][1];
        OpM13Vec[i] = DummyRAVar_operator[i][0][2];
        OpM21Vec[i] = DummyRAVar_operator[i][1][0];
        OpM22Vec[i] = DummyRAVar_operator[i][1][1];
        OpM23Vec[i] = DummyRAVar_operator[i][1][2];
        OpM31Vec[i] = DummyRAVar_operator[i][2][0];
        OpM32Vec[i] = DummyRAVar_operator[i][2][1];
        OpM33Vec[i] = DummyRAVar_operator[i][2][2];
    }
}

class CAECyldrField
{
private:
    int m_NoOfCyldr;                                // Number of Cylinders
    int m_NoOfFieldPoint;                           // Number of field data points
    double m_Mu;                                    // permiability
    int m_Phi_Steps;                                // number of steps for angle phi
    double m_Acc_r;                                 // accuracy for r

    // Ir2
    double Solve_Ir2(double Ro, double FP_R, double Phi1, double Phi2, double w1, double w2)
    {
        // not equal to zero
        auto lmbda_I_r2_11_n0 = [FP_R, Phi1, w1](double rN)
        {
            return -(rN/FP_R)* log(w1 + sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w1, 2) - 2* FP_R* rN* cos(Phi1)) );
        };

        auto lmbda_I_r2_12_n0 = [FP_R, Phi1, w2](double rN)
        {
            return -(rN/FP_R)* log(w2 + sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w2, 2) - 2* FP_R* rN* cos(Phi1)) );
        };

        auto lmbda_I_r2_21_n0 = [FP_R, Phi2, w1](double rN)
        {
            return -(rN/FP_R)* log(w1 + sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w1, 2) - 2* FP_R* rN* cos(Phi2)) );
        };

        auto lmbda_I_r2_22_n0 = [FP_R, Phi2, w2](double rN)
        {
            return -(rN/FP_R)* log(w2 + sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w2, 2) - 2* FP_R* rN* cos(Phi2)) );
        };

        // equal to zero
        auto lmbda_I_r2_0 = [Ro](double Phi_j, double wk)
        {
            return - cos(Phi_j)* wk* asinh(Ro/ wk);
        };

        int StepsR = 1 + (Ro - 0)/ m_Acc_r;

        double I_r2;

        if (FP_R == 0)
        {
            I_r2 =   lmbda_I_r2_0(Phi1, w1)
                   - lmbda_I_r2_0(Phi1, w2)
                   - lmbda_I_r2_0(Phi2, w1)
                   + lmbda_I_r2_0(Phi2, w2);
        }
        else
        {
            I_r2 =   Trapz(lmbda_I_r2_11_n0, StepsR, 0, Ro)
                   - Trapz(lmbda_I_r2_12_n0, StepsR, 0, Ro)
                   - Trapz(lmbda_I_r2_21_n0, StepsR, 0, Ro)
                   + Trapz(lmbda_I_r2_22_n0, StepsR, 0, Ro);
        }
        return I_r2;
    }

    // Ir13
    double Solve_Ir13(double Ro, double FP_R, double Phi1, double Phi2, double w1, double w2)
    {
        auto lmbda_Ir13_1 = [FP_R, w1](double rN, double PhiN)
        {
            return w1* (FP_R* rN - pow(rN, 2)* cos(PhiN) ) / 
                    ( ( pow(FP_R, 2) + pow(rN, 2) - 2* FP_R* rN* cos(PhiN) )* 
                        sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w1, 2) - 2* FP_R* rN* cos(PhiN)) );
        };

        auto lmbda_Ir13_2 = [FP_R, w2](double rN, double PhiN)
        {
            return w2* (FP_R* rN - pow(rN, 2)* cos(PhiN) ) / 
                    ( ( pow(FP_R, 2) + pow(rN, 2) - 2* FP_R* rN* cos(PhiN) )* 
                        sqrt(pow(FP_R, 2) + pow(rN, 2) + pow(w2, 2) - 2* FP_R* rN* cos(PhiN)) );
        };

        int StepsR = 1 + (Ro - 0)/ m_Acc_r;

        double I_r13 = - Trapz2D(lmbda_Ir13_1, StepsR, 0, Ro, m_Phi_Steps, Phi1, Phi2)
                       + Trapz2D(lmbda_Ir13_2, StepsR, 0, Ro, m_Phi_Steps, Phi1, Phi2);

        return I_r13;
    }

    // Magnetic field calculation for one data point due to one filament only. One to one map.
    tuple<double, double, double> CyldrMagnetic(double X,     double Y,     double Z,                // coordinates of field point
                                                double Ro,    double L,     double J,                // radius, length, and current density
                                                double Xo,    double Yo,    double Zo,               // center
                                                double OpM11, double OpM12, double OpM13,            // operator variable
                                                double OpM21, double OpM22, double OpM23,
                                                double OpM31, double OpM32, double OpM33)
    {
        // Inverse transformation
        double FP_X = OpM11* (X-Xo) + OpM21* (Y-Yo) + OpM31* (Z-Zo);
        double FP_Y = OpM12* (X-Xo) + OpM22* (Y-Yo) + OpM32* (Z-Zo);
        double FP_Z = OpM13* (X-Xo) + OpM23* (Y-Yo) + OpM33* (Z-Zo);

        double FP_R, FP_Phy;                                                      // FP means Field Point
        tie(FP_R, FP_Phy) = Cart2Pol(FP_X, FP_Y);                                 // field point in polar fashion

        // constants
        double Phi1  = - Pi - FP_Phy;
        double Phi2  =   Pi - FP_Phy;
        double w1    = - 0.5* L - FP_Z;
        double w2    =   0.5* L - FP_Z;

        double I_r2  = Solve_Ir2(Ro, FP_R, Phi1, Phi2, w1, w2);
        double I_r13 = Solve_Ir13(Ro, FP_R, Phi1, Phi2, w1, w2);

        double Hx    = cos(FP_Phy)* I_r2 - sin(FP_Phy)* I_r13;
        double Hy    = sin(FP_Phy)* I_r2 + cos(FP_Phy)* I_r13;
        double Hz    = 0;

        double Bix   = (m_Mu* J/ (4* Pi))* Hx;
        double Biy   = (m_Mu* J/ (4* Pi))* Hy;
        double Biz   = (m_Mu* J/ (4* Pi))* Hz;

        // Reverse transformation
        double Bx = OpM11* Bix + OpM12* Biy + OpM13* Biz;
        double By = OpM21* Bix + OpM22* Biy + OpM23* Biz;
        double Bz = OpM31* Bix + OpM32* Biy + OpM33* Biz;
        
        return make_tuple(Bx, By, Bz);
    }

    // Multiple filaments and one data point. Many to one map.
    tuple<double, double, double> ManyCyldr(double X,         double Y,         double Z,
                                            double *RoVec,    double *LVec,     double *JVec,
                                            double *XoVec,    double *YoVec,    double *ZoVec,
                                            double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                                            double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                                            double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        double Bx = 0.0, By = 0.0, Bz = 0.0;    // Field due to multiple filaments
        double bx, by, bz;                      // Dummy variables

        for(int i = 0; i < m_NoOfCyldr; i++)
        {
            tie(bx, by, bz) = CyldrMagnetic(X, Y, Z,                                  // remember we need only one data point
                                            RoVec[i],    LVec[i],     JVec[i],
                                            XoVec[i],    YoVec[i],    ZoVec[i],
                                            OpM11Vec[i], OpM12Vec[i], OpM13Vec[i],
                                            OpM21Vec[i], OpM22Vec[i], OpM23Vec[i],
                                            OpM31Vec[i], OpM32Vec[i], OpM33Vec[i]);

            Bx += bx;
            By += by;
            Bz += bz;
        }

        return make_tuple(Bx, By, Bz);
    }

public:

    CAECyldrField(int noOfCyldr, int noOfFieldPoint, double mu, int phi_Steps, double acc_r)
    {
        m_NoOfCyldr      = noOfCyldr;
        m_NoOfFieldPoint = noOfFieldPoint;
        m_Mu             = mu;
        m_Phi_Steps      = phi_Steps;
        m_Acc_r          = acc_r;
    }

    // Magnetic field at each data point
    void CyldrMagneticField(double *BxAll,    double *ByAll,    double *BzAll,
                            double *X,        double *Y,        double *Z,
                            double *RoVec,    double *LVec,     double *JVec,
                            double *XoVec,    double *YoVec,    double *ZoVec,
                            double *OpM11Vec, double *OpM12Vec, double *OpM13Vec,
                            double *OpM21Vec, double *OpM22Vec, double *OpM23Vec,
                            double *OpM31Vec, double *OpM32Vec, double *OpM33Vec)
    {
        cout << "Calculating magnetic field due to cylinder." << endl;

        for(int i = 0; i < m_NoOfFieldPoint; i++)
        {
            tie(BxAll[i], ByAll[i], BzAll[i]) = ManyCyldr(X[i], Y[i], Z[i],
                                                          RoVec,    LVec,     JVec,
                                                          XoVec,    YoVec,    ZoVec,
                                                          OpM11Vec, OpM12Vec, OpM13Vec,
                                                          OpM21Vec, OpM22Vec, OpM23Vec,
                                                          OpM31Vec, OpM32Vec, OpM33Vec);

            cout << i + 1 << "   ";
            if( (i+1)%10 == 0 ){cout << "\n";}
        }
        cout << "\n";
    }

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int main()
{
    double mu;
    int NoOfFieldPoints, NoOfFilament, NoOfSlab, NoOfArc, NoOfRectArc, NoOfCylndr, NoOfHelix;
    
    fstream file("Bridge/ImportantVar.txt", ios_base::in);
    file >> mu >> NoOfFieldPoints >> NoOfFilament >> NoOfSlab >> NoOfArc >> NoOfRectArc >> NoOfCylndr >> NoOfHelix;
    file.close();

    // Total Magnetic field array
    double TotalBx[NoOfFieldPoints], TotalBy[NoOfFieldPoints], TotalBz[NoOfFieldPoints];
    InitVec(NoOfFieldPoints, TotalBx, TotalBy, TotalBz); // Initialize them

    // Load field point
    double FieldXvec[NoOfFieldPoints], FieldYvec[NoOfFieldPoints], FieldZvec[NoOfFieldPoints];
    LoadFieldPoints(NoOfFieldPoints, FieldXvec, FieldYvec, FieldZvec);

    // Points in Polar
    double FieldRvec[NoOfFieldPoints], FieldThetavec[NoOfFieldPoints];
    for(int i = 0; i < NoOfFieldPoints; i++)
    {
        tie(FieldRvec[i], FieldThetavec[i]) = Cart2Pol(FieldXvec[i], FieldYvec[i]);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfFilament > 0)
    {
        double LengthVec[NoOfFilament], IVec[NoOfFilament];                                 // lenght and current in each filament. Data is in vector formant
        double XoVec[NoOfFilament],     YoVec[NoOfFilament],    ZoVec[NoOfFilament];        // center coordinate of each filament
        double OpM11Vec[NoOfFilament],  OpM12Vec[NoOfFilament], OpM13Vec[NoOfFilament],     // operator elements
               OpM21Vec[NoOfFilament],  OpM22Vec[NoOfFilament], OpM23Vec[NoOfFilament], 
               OpM31Vec[NoOfFilament],  OpM32Vec[NoOfFilament], OpM33Vec[NoOfFilament];
        
        // Load data
        LoadFilamentData(NoOfFilament,
                         LengthVec, IVec,
                         XoVec,     YoVec,    ZoVec,
                         OpM11Vec,  OpM12Vec, OpM13Vec,
                         OpM21Vec,  OpM22Vec, OpM23Vec,
                         OpM31Vec,  OpM32Vec, OpM33Vec);

        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                  // Initialize them

        CAEFilamentField CalculateFilament = CAEFilamentField(NoOfFilament, NoOfFieldPoints, mu);
        
        CalculateFilament.FilamentMagneticField(dumBx,     dumBy,     dumBz,
                                                FieldXvec, FieldYvec, FieldZvec,
                                                LengthVec, IVec, 
                                                XoVec,     YoVec,     ZoVec,
                                                OpM11Vec,  OpM12Vec,  OpM13Vec,
                                                OpM21Vec,  OpM22Vec,  OpM23Vec,
                                                OpM31Vec,  OpM32Vec,  OpM33Vec);

        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfSlab > 0)
    {
        double JVec[NoOfSlab];                                                  // Current density in each slab
        double XoVec[NoOfSlab],    YoVec[NoOfSlab],    ZoVec[NoOfSlab];         // center coordinate of each slab
        double XVerOrg0[NoOfSlab], YVerOrg0[NoOfSlab], ZVerOrg0[NoOfSlab],      // two extreme points of each slab
               XVerOrg6[NoOfSlab], YVerOrg6[NoOfSlab], ZVerOrg6[NoOfSlab];
        double OpM11Vec[NoOfSlab], OpM12Vec[NoOfSlab], OpM13Vec[NoOfSlab],      // operator for each slab
               OpM21Vec[NoOfSlab], OpM22Vec[NoOfSlab], OpM23Vec[NoOfSlab], 
               OpM31Vec[NoOfSlab], OpM32Vec[NoOfSlab], OpM33Vec[NoOfSlab];

        LoadSlabData(NoOfSlab, JVec,
                     XVerOrg0, YVerOrg0, ZVerOrg0,
                     XVerOrg6, YVerOrg6, ZVerOrg6,
                     XoVec,    YoVec,    ZoVec,
                     OpM11Vec, OpM12Vec, OpM13Vec,
                     OpM21Vec, OpM22Vec, OpM23Vec,
                     OpM31Vec, OpM32Vec, OpM33Vec);

        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                  // Initialize them
    
        CAESlabField CalculateSlab = CAESlabField(NoOfSlab, NoOfFieldPoints, mu);
        
        CalculateSlab.SlabMagneticField(dumBx,     dumBy,     dumBz,
                                        FieldXvec, FieldYvec, FieldZvec,
                                        XVerOrg0,  YVerOrg0,  ZVerOrg0,
                                        XVerOrg6,  YVerOrg6,  ZVerOrg6,
                                        JVec, 
                                        XoVec,     YoVec,     ZoVec,
                                        OpM11Vec,  OpM12Vec,  OpM13Vec,
                                        OpM21Vec,  OpM22Vec,  OpM23Vec,
                                        OpM31Vec,  OpM32Vec,  OpM33Vec);

        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfArc > 0)
    {
        double IVec[NoOfArc],     a0Vec[NoOfArc];                            // Current in each arc
        double XoVec[NoOfArc],    YoVec[NoOfArc],    ZoVec[NoOfArc];         // center coordinate of each arc
        double Phy1Vec[NoOfArc],  Phy2Vec[NoOfArc];
        double OpM11Vec[NoOfArc], OpM12Vec[NoOfArc], OpM13Vec[NoOfArc],      // operator for each arc
               OpM21Vec[NoOfArc], OpM22Vec[NoOfArc], OpM23Vec[NoOfArc], 
               OpM31Vec[NoOfArc], OpM32Vec[NoOfArc], OpM33Vec[NoOfArc];
        
        // Load data
        LoadArcData(NoOfArc, 
                    IVec,     a0Vec,
                    Phy1Vec,  Phy2Vec,
                    XoVec,    YoVec,    ZoVec,
                    OpM11Vec, OpM12Vec, OpM13Vec,
                    OpM21Vec, OpM22Vec, OpM23Vec,
                    OpM31Vec, OpM32Vec, OpM33Vec);

        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                  // Initialize them

        CAEArcField CalculateArc = CAEArcField(NoOfArc, NoOfFieldPoints, mu);

        CalculateArc.ArcMagneticField(dumBx,     dumBy,     dumBz,
                                      FieldXvec, FieldYvec, FieldZvec,
                                      a0Vec,     IVec,
                                      Phy1Vec,   Phy2Vec,
                                      XoVec,     YoVec,     ZoVec,
                                      OpM11Vec,  OpM12Vec,  OpM13Vec,
                                      OpM21Vec,  OpM22Vec,  OpM23Vec,
                                      OpM31Vec,  OpM32Vec,  OpM33Vec);

        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfRectArc > 0)
    {
        double RInVec[NoOfRectArc],   ROutVec[NoOfRectArc];                              // Inner and Outer radius of each rectangular arc
        double TcknsVec[NoOfRectArc], JVec[NoOfRectArc];                                 // Thickness and current density in rectangular arc
        double Phy1Vec[NoOfRectArc],  Phy2Vec[NoOfRectArc];                              // Initial and Final angle of each rectangular arc
        double XoVec[NoOfRectArc],    YoVec[NoOfRectArc],    ZoVec[NoOfRectArc];         // center coordinate of each ractangular arc
        double OpM11Vec[NoOfRectArc], OpM12Vec[NoOfRectArc], OpM13Vec[NoOfRectArc],      // operator elements
               OpM21Vec[NoOfRectArc], OpM22Vec[NoOfRectArc], OpM23Vec[NoOfRectArc], 
               OpM31Vec[NoOfRectArc], OpM32Vec[NoOfRectArc], OpM33Vec[NoOfRectArc];

        // Load data
        LoadRectArcData(NoOfRectArc, 
                        RInVec,   ROutVec,
                        TcknsVec, JVec,
                        Phy1Vec,  Phy2Vec,
                        XoVec,    YoVec,    ZoVec,
                        OpM11Vec, OpM12Vec, OpM13Vec,
                        OpM21Vec, OpM22Vec, OpM23Vec,
                        OpM31Vec, OpM32Vec, OpM33Vec);

        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                   // Initialize them


        // For the integrations with respect to Phi (relative angle), 
        // the programmer has chosen "number of steps" method in numerical integration
        // because Phi max to max varies from 0 to 2*Pi, so 500-1000 steps are egnough.
        //
        // While in the case of radial integration, the programmer has chosen "accuracy" method
        // as r can vary from 0 to infinity. So accuracy is the vice option.

        CAERectArcField CalculateRectArc = CAERectArcField(NoOfRectArc, NoOfFieldPoints, mu, 500, 0.005);

        CalculateRectArc.RectArcMagneticField(dumBx,     dumBy,     dumBz,
                                              FieldXvec, FieldYvec, FieldZvec,
                                              RInVec,    ROutVec,
                                              TcknsVec,  JVec,
                                              Phy1Vec,   Phy2Vec,
                                              XoVec,     YoVec,     ZoVec,
                                              OpM11Vec,  OpM12Vec,  OpM13Vec,
                                              OpM21Vec,  OpM22Vec,  OpM23Vec,
                                              OpM31Vec,  OpM32Vec,  OpM33Vec);


        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfCylndr > 0)
    {
        double RoVec[NoOfCylndr],    LVec[NoOfCylndr],     JVec[NoOfCylndr];          // radius and length of cylinder
        double XoVec[NoOfCylndr],    YoVec[NoOfCylndr],    ZoVec[NoOfCylndr];         // center coordinate of each ractangular arc
        double OpM11Vec[NoOfCylndr], OpM12Vec[NoOfCylndr], OpM13Vec[NoOfCylndr],      // operator elements
               OpM21Vec[NoOfCylndr], OpM22Vec[NoOfCylndr], OpM23Vec[NoOfCylndr], 
               OpM31Vec[NoOfCylndr], OpM32Vec[NoOfCylndr], OpM33Vec[NoOfCylndr];

        LoadCyldrData(NoOfCylndr, 
                      RoVec,    LVec,     JVec,
                      XoVec,    YoVec,    ZoVec,
                      OpM11Vec, OpM12Vec, OpM13Vec,
                      OpM21Vec, OpM22Vec, OpM23Vec,
                      OpM31Vec, OpM32Vec, OpM33Vec);

        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                  // Initialize them
    
        CAECyldrField CalculateCyldr = CAECyldrField(NoOfCylndr, NoOfFieldPoints, mu, 500, 0.005);

        CalculateCyldr.CyldrMagneticField(dumBx,     dumBy,     dumBz,
                                          FieldXvec, FieldYvec, FieldZvec,
                                          RoVec,     LVec,      JVec,
                                          XoVec,     YoVec,     ZoVec,
                                          OpM11Vec,  OpM12Vec,  OpM13Vec,
                                          OpM21Vec,  OpM22Vec,  OpM23Vec,
                                          OpM31Vec,  OpM32Vec,  OpM33Vec);

        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (NoOfHelix > 0)
    {


        // Dummy magnetic field
        double dumBx[NoOfFieldPoints], dumBy[NoOfFieldPoints], dumBz[NoOfFieldPoints];
        InitVec(NoOfFieldPoints, dumBx, dumBy, dumBz);                                  // Initialize them
    
        // Magnetic field calculation


        // Add this magnetic field to the global magnetic field
        AddVecs(NoOfFieldPoints, TotalBx, dumBx);
        AddVecs(NoOfFieldPoints, TotalBy, dumBy);
        AddVecs(NoOfFieldPoints, TotalBz, dumBz);
    }


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    WriteDataToCSV(NoOfFieldPoints, FieldXvec, FieldYvec, FieldZvec, FieldRvec, FieldThetavec, TotalBx, TotalBy, TotalBz);

}
