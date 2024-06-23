#include <math.h>
#include <mex.h>

double rem(double x, double y);

// Global variables
double pi = 3.1415926535897931;
double rad = pi/180;
int slsr[60][6] = {
					{         0,         0,         1,         0,   6288774, -20905355},
					{         2,         0,        -1,         0,   1274027,  -3699111},
					{         2,         0,         0,         0,    658314,  -2955968},
					{         0,         0,         2,         0,    213618,   -569925},
					{         0,         1,         0,         0,   -185116,     48888},
					{         0,         0,         0,         2,   -114332,     -3149},
					{         2,         0,        -2,         0,     58793,    246158},
					{         2,        -1,        -1,         0,     57066,   -152138},
					{         2,         0,         1,         0,     53322,   -170733},
					{         2,        -1,         0,         0,     45758,   -204586},
					{         0,         1,        -1,         0,    -40923,   -129620},
					{         1,         0,         0,         0,    -34720,    108743},
					{         0,         1,         1,         0,    -30383,    104755},
					{         2,         0,         0,        -2,     15327,     10321},
					{         0,         0,         1,         2,    -12528,         0},
					{         0,         0,         1,        -2,     10980,     79661},
					{         4,         0,        -1,         0,     10675,    -34782},
					{         0,         0,         3,         0,     10034,    -23210},
					{         4,         0,        -2,         0,      8548,    -21636},
					{         2,         1,        -1,         0,     -7888,     24208},
					{         2,         1,         0,         0,     -6766,     30824},
					{         1,         0,        -1,         0,     -5163,     -8379},
					{         1,         1,         0,         0,      4987,    -16675},
					{         2,        -1,         1,         0,      4036,    -12831},
					{         2,         0,         2,         0,      3994,    -10445},
					{         4,         0,         0,         0,      3861,    -11650},
					{         2,         0,        -3,         0,      3665,     14403},
					{         0,         1,        -2,         0,     -2689,     -7003},
					{         2,         0,        -1,         2,     -2602,         0},
					{         2,        -1,        -2,         0,      2390,     10056},
					{         1,         0,         1,         0,     -2348,      6322},
					{         2,        -2,         0,         0,      2236,     -9884},
					{         0,         1,         2,         0,     -2120,      5751},
					{         0,         2,         0,         0,     -2069,         0},
					{         2,        -2,        -1,         0,      2048,     -4950},
					{         2,         0,         1,        -2,     -1773,      4130},
					{         2,         0,         0,         2,     -1595,         0},
					{         4,        -1,        -1,         0,      1215,     -3958},
					{         0,         0,         2,         2,     -1110,         0},
					{         3,         0,        -1,         0,      -892,      3258},
					{         2,         1,         1,         0,      -810,      2616},
					{         4,        -1,        -2,         0,       759,     -1897},
					{         0,         2,        -1,         0,      -713,     -2117},
					{         2,         2,        -1,         0,      -700,      2354},
					{         2,         1,        -2,         0,       691,         0},
					{         2,        -1,         0,        -2,       596,         0},
					{         4,         0,         1,         0,       549,     -1423},
					{         0,         0,         4,         0,       537,     -1117},
					{         4,        -1,         0,         0,       520,     -1571},
					{         1,         0,        -2,         0,      -487,     -1739},
					{         2,         1,         0,        -2,      -399,         0},
					{         0,         0,         2,        -2,      -381,     -4421},
					{         1,         1,         1,         0,       351,         0},
					{         3,         0,        -2,         0,      -340,         0},
					{         4,         0,        -3,         0,       330,         0},
					{         2,        -1,         2,         0,       327,         0},
					{         0,         2,         1,         0,      -323,      1165},
					{         1,         1,        -1,         0,       299,         0},
					{         2,         0,         3,         0,       294,         0},
					{         2,         0,        -1,        -2,         0,      8752}
					};
int sb[60][5] =	{
					{         0,         0,         0,         1,   5128122},
					{         0,         0,         1,         1,    280602},
					{         0,         0,         1,        -1,    277693},
					{         2,         0,         0,        -1,    173237},
					{         2,         0,        -1,         1,     55413},
					{         2,         0,        -1,        -1,     46271},
					{         2,         0,         0,         1,     32573},
					{         0,         0,         2,         1,     17198},
					{         2,         0,         1,        -1,      9266},
					{         0,         0,         2,        -1,      8822},
					{         2,        -1,         0,        -1,      8216},
					{         2,         0,        -2,        -1,      4324},
					{         2,         0,         1,         1,      4200},
					{         2,         1,         0,        -1,     -3359},
					{         2,        -1,        -1,         1,      2463},
					{         2,        -1,         0,         1,      2211},
					{         2,        -1,        -1,        -1,      2065},
					{         0,         1,        -1,        -1,     -1870},
					{         4,         0,        -1,        -1,      1828},
					{         0,         1,         0,         1,     -1794},
					{         0,         0,         0,         3,     -1749},
					{         0,         1,        -1,         1,     -1565},
					{         1,         0,         0,         1,     -1491},
					{         0,         1,         1,         1,     -1475},
					{         0,         1,         1,        -1,     -1410},
					{         0,         1,         0,        -1,     -1344},
					{         1,         0,         0,        -1,     -1335},
					{         0,         0,         3,         1,      1107},
					{         4,         0,         0,        -1,      1021},
					{         4,         0,        -1,         1,       833},
					{         0,         0,         1,        -3,       777},
					{         4,         0,        -2,         1,       671},
					{         2,         0,         0,        -3,       607},
					{         2,         0,         2,        -1,       596},
					{         2,        -1,         1,        -1,       491},
					{         2,         0,        -2,         1,      -451},
					{         0,         0,         3,        -1,       439},
					{         2,         0,         2,         1,       422},
					{         2,         0,        -3,        -1,       421},
					{         2,         1,        -1,         1,      -366},
					{         2,         1,         0,         1,      -351},
					{         4,         0,         0,         1,       331},
					{         2,        -1,         1,         1,       315},
					{         2,        -2,         0,        -1,       302},
					{         0,         0,         1,         3,      -283},
					{         2,         1,         1,        -1,      -229},
					{         1,         1,         0,        -1,       223},
					{         1,         1,         0,         1,       223},
					{         0,         1,        -2,        -1,      -220},
					{         2,         1,        -1,        -1,      -220},
					{         1,         0,         1,         1,      -185},
					{         2,        -1,        -2,        -1,       181},
					{         0,         1,         2,         1,      -177},
					{         4,         0,        -2,        -1,       176},
					{         4,        -1,        -1,        -1,       166},
					{         1,         0,         1,        -1,      -164},
					{         4,         0,         1,        -1,       132},
					{         1,         0,        -1,        -1,      -119},
					{         4,        -1,         0,        -1,       115},
					{         2,        -2,         0,         1,       107}
					};


void lunarpos(double aapp[],double dapp[],double r[],double lam[],double beta[],const double t[],const double eps0[],const double deleps[],const double delpsi[],int row,int col)
{
	if(col > row)
	{
		row = col;
	}

	for(int i = 0; i < row; i++)
	{
		// Local variables
		double Ls = 0;
		double D  = 0;
		double M  = 0;
		double Ms = 0;
		double F  = 0;
		double A1 = 0;
		double A2 = 0;
		double A3 = 0;

		double E  = 0;
		double Eq = 0;

		double scarg1 = 0;
		double scarg2 = 0;
		double scoef1 = 0;
		double ccoef1 = 0;
		double scoef2 = 0;
		double sigl   = 0;
		double sigr   = 0;
		double sigb   = 0;

		double lams = 0;
		double eps  = 0;
		double b2   = 0;

		//-------------------------------------------------------
		// the angles Ls,D,M,Ms,F from Meeus S.336/337, all in [rad]
		//-------------------------------------------------------
		Ls = rad*(-pow(t[i],4)/65194000  + pow(t[i],3)/538841   + -pow(t[i],2)*0.0013268 + t[i]*481267.88134236 + 218.3164591);
		D  = rad*(-pow(t[i],4)/113065000 + pow(t[i],3)/545868   + -pow(t[i],2)*0.0016300 + t[i]*445267.1115168  + 297.8502042);
		M  = rad*( 0                     + pow(t[i],3)/24490000 + -pow(t[i],2)*0.0001536 + t[i]*35999.0502909   + 357.5291092);
		Ms = rad*(-pow(t[i],4)/14712000  + pow(t[i],3)/69699    + -pow(t[i],2)*0.0089970 + t[i]*477198.8676313  + 134.9634114);
		F  = rad*( pow(t[i],4)/863310000 + -pow(t[i],3)/3526000 + -pow(t[i],2)*0.0034029 + t[i]*483202.0175273  + 93.2720993);

		A1 = rad*(119.75 +    131.849*t[i]);
		A2 = rad*( 53.09 + 479264.290*t[i]);
		A3 = rad*(313.45 + 481266.484*t[i]);


		//----------------------------------------------------
		// computation of the sums sigl, sigr, sigb
		// slsr.mat and sb.mat contain tables from Meeus
		//----------------------------------------------------
		// periodic terms for long. and dist. to Moon
		// periodic terms for latitude of Moon
		E  = 1 - 0.002516*t[i] - 0.0000074*pow(t[i],2);
		Eq = pow(E,2);

		sigl = 0;
		sigr = 0;
		sigb = 0;
		for(int k = 0; k < 60; k++)
		{

			scarg1 = slsr[k][0]*D + slsr[k][1]*M + slsr[k][2]*Ms + slsr[k][3]*F;
			scarg2 = sb[k][0]*D + sb[k][1]*M + sb[k][2]*Ms + sb[k][3]*F;

			if(abs(slsr[k][1]) == 1)
			{
				scoef1 = slsr[k][4]*E;
				ccoef1 = slsr[k][5]*E;
				scoef2 = sb[k][4]*E;
			}
			else if(abs(slsr[k][1]) == 2)
			{
				scoef1 = slsr[k][4]*Eq;
				ccoef1 = slsr[k][5]*Eq;
				scoef2 = sb[k][4]*Eq;
			}
			else
			{
				scoef1 = slsr[k][4];
				ccoef1 = slsr[k][5];
				scoef2 = sb[k][4];
			}

			sigl += scoef1*sin(scarg1);
			sigr += ccoef1*cos(scarg1);
			sigb += scoef2*sin(scarg2);
		}
		sigl = (sigl + 3958*sin(A1) + 1962*sin(Ls-F) + 318*sin(A2))/1000000;
		sigr = sigr/1000;
		sigb = (sigb - 2235*sin(Ls) + 382*sin(A3) + 175*sin(A1-F) + 175*sin(A1+F) + 127*sin(Ls-Ms) - 115*sin(Ls+Ms))/1000000;

		//-------------------------
		// coordinates lam,beta,r
		//-------------------------
		lam[i]  = rem(rem(Ls*180/pi+sigl,360)+360,360);
		beta[i] = sigb;
		r[i]    = (385000.56 + sigr)*1000;

		//--------------------------------------------
		// apparent right ascension and declination
		//--------------------------------------------
		// nutation in longitude (accuracy ~ 1")
		// in Meeus more exact formulae (e.g. p.339)
		lams = (lam[i]+delpsi[i])*pi/180;
		eps  = (eps0[i]+deleps[i])*pi/180;
		b2   = beta[i]*pi/180;
		aapp[i] = atan2((sin(lams)*cos(eps)-tan(b2)*sin(eps)),cos(lams));
		aapp[i] = rem(aapp[i]*180/pi+360,360);
		dapp[i] = asin(sin(b2)*cos(eps)+cos(b2)*sin(eps)*sin(lams));
		dapp[i] = dapp[i]*180/pi;

	}
}

double rem(double x, double y)
{
	return x - floor(x/y)*y;
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	// Inputs
	double *t;
	double *eps0;
	double *deleps;
	double *delpsi;

	// Outputs
	double *aapp;
	double *dapp;
	double *r;
	double *lam;
	double *beta;

	// Get number of rows
	int row = mxGetM(prhs[0]);
	int col = mxGetN(prhs[0]);

	// Create a 1-by-N matrix for the return argument
	plhs[0] = mxCreateDoubleMatrix(row,col,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(row,col,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(row,col,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(row,col,mxREAL);
	plhs[4] = mxCreateDoubleMatrix(row,col,mxREAL);

	// Assign a pointer to the intput
	t      = mxGetPr(prhs[0]);
	eps0   = mxGetPr(prhs[1]);
	deleps = mxGetPr(prhs[2]);
	delpsi = mxGetPr(prhs[3]);

	// Assign a pointer to the output
	aapp = mxGetPr(plhs[0]);
	dapp = mxGetPr(plhs[1]);
	r    = mxGetPr(plhs[2]);
	lam  = mxGetPr(plhs[3]);
	beta = mxGetPr(plhs[4]);

	// Call the lunarpos function
	lunarpos(aapp,dapp,r,lam,beta,t,eps0,deleps,delpsi,row,col);
}