#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <math.h>
#include <complex>
#include <time.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <assert.h>

/*constants*/
#define EPS_0 8.85418782e-12  	// F/m, vacuum permittivity
#define K	1.38065e-23			// J/K, Boltzmann constant
#define ME 9.10938215e-31		// kg, electron mass
#define QE 1.602176565e-19		// C, elementary charge
#define AMU  1.660538921e-27	// kg, atomic mass unit
#define EV_TO_K	11604.52		// 1eV in Kelvin, QE/K
#define PI 3.14159265//

#define quiet 0
#define MEM_FAC 4

/* Data structure for particle storage **/
struct Particle
{
	double v[3];			/*velocity*/
	double r[3];			/*position for diffusion calculation*/
	double mass;			/*particle mass in kg*/
	double charge;			/*particle charge in Coulomb*/
	int    type;            /* Species type*/
};


struct ParticleList
{
	std::vector<unsigned long> particleIdentifier;
};

/* Data structure to hold species information*/
struct Particles
{

	unsigned long np;					    /*number of particles*/
	unsigned long np_alloc;			    /*size of the allocated data array*/
	int numTypes;				/*number of types*/
	Particle *part;			    /*array holding particles*/
	ParticleList *typeList; /*list of particles of a particular type*/

	/*Functions*/
	double Vrel(unsigned long partid1, unsigned long partid2); //Function to calculate Vrel.
	double Erel(unsigned long partid1, unsigned long partid2);
	double Etot(unsigned long partid1, unsigned long partid2);
	int RandOfType(int type, unsigned long *pid);
};

struct CrossSec
{

	std::vector<double> sigma;			/*array holding particles*/
	std::vector<double> energy;		/*array holding particles*/
	int ParticleType1;      /*Identifier for particle type in collision*/
	int ParticleType2;      /*Identifier for particle type in collision*/
	int ProductType1;
	int ProductType2;
	unsigned long count; /*number of interactions of this type that have*/
	std::string collisionType;
	double DeltaE;
	/*Functions*/
	double SigmaOfE(double PartEnergy, bool extrapolate);

};

struct CrossSections
{
	int NCrossSec;
	CrossSec *csec;
	double *VSigmaMax; /*Store samples of max(v sigma)*/
    double *remainders; /*Store Npair remainders*/
};

struct DifferentialCrossSection {
    std::string modelName;
    std::vector<double> params;
};

struct InputData
{
	std::vector<std::string> collisionFilename;
	std::vector<int> particle1, particle2,product1,product2;
	std::vector<std::string> collisionType;
	std::vector<double> DeltaE;
    std::vector<DifferentialCrossSection> DCS;

    std::vector<std::string> velocityLoadFilenames;

	double E;   //Electric Field
	double Efreq; //Frequency of E field oscillation
	double Emag; //E field Magnitude
	double B[3];   // Magnetic Field
	int NUM_TS; //Number of time steps
	double DT;				// time step size
	double BOX_LENGTH;
	int NSigmaMaxSamples;
	int NUM_TYPES;
	std::vector<unsigned long> Nparticles;      
	double *Tparticles;  //temperature in eV
	double *Vparticles;  //flow shift in eV
	double *Mparticles;  //mass in amu
	double *Qparticles;  //charge in |e|
    std::string CollOrder;
    double MEM;          // Memory allocation request
	
	//output
    std::vector<std::vector<int>> DumpVDFParams;
    int Slurm; // Slurm exit flag
    int OS; // Output stride
    int CR; // Carry remainder to next time step in NTC Npairs eq.
    int Extrapolate; // Extrapolate cross sections
};

/** FUNCTION PROTOTYPES **/
double rnd();
void ReportStatus(Particles *AllParticle, CrossSections *AllCrossSections, InputData *Parameters, int step);
double SampleVel(double T, double mass);

// Include in future particle class
void AddParticle(InputData *SimulationParam, Particles *particle_list,
                 double rx, double ry, double rz, double vx, double vy, double vz,
                 double charge, double mass, int type, int step);
void InitializeParticles(Particles *particle_list, InputData *SimulationParam);
void Accelerate(Particles *particle_list, InputData *SimulationParam);
void AccelerateBoris(Particles *particle_list, InputData *SimulationParam);
void UpdateDisplacement(Particles *particle_list,InputData *SimulationParam, int type);
double CalculateMeanDisplacement(Particles *particle_list,int axis, int type);
double ComputeKE(Particles *particle_list, int type);	/*computes kinetic energy*/
double ComputeMeanE(Particles *particle_list, InputData *SimulationParam,int type);	/*computes mean energy*/
double ComputeVelocityMoment(Particles *particle_list, int type, int axis);	
double ComputeTemperature(Particles *particle_list, InputData *SimulationParam,int type,int axis);
void DumpParticles(Particles *particle_list, InputData *Parameters, int step);
int SlurmExit(int i, int NS);

// Include in future cross section class
void ReadCrossSection(InputData *SimulationParam, CrossSections *CrossSecList);
void SampleVSigmaMax(Particles *particle_list, CrossSections *CrossSecList,int samples);
void DumpReactionCounts(CrossSections *CrossSecList,int step);

void DSMCcollideNTC(Particles *particle_list, CrossSections *CrossSecList, InputData *SimulationParam, int step);

// Collision models
void ElasticCollision(Particles *particle_list, unsigned long RandomParticleType1,
                      unsigned long RandomParticleType2, bool FixedHeavyParticle,
                      DifferentialCrossSection dcs);
void InelasticCollision(Particles *particle_list, unsigned long RandomParticleType1,
                        unsigned long RandomParticleType2, double DeltaE, bool FixedHeavyParticle);
void CXCollision(Particles *particle_list, unsigned long RandomParticleType1, unsigned long RandomParticleType2);
// Helpers
unsigned int ReadTime(std::string g);
void RotateFrame(double *v, double ct, double st, double cp, double sp);
void Project(double *v, double cos_chi, double sin_chi, double cos_eta, double sin_eta);
void vadd_scale(double *v, double *v1, double *v2, double scale);

// Simulation parameters
void ReadInput(std::string inputFileName, InputData *SimulationParam);
void UpdateEField(InputData *SimulationParam, int step);

// Differential distribution models
double ParkElasticCollisionCosX(double eps, std::vector<double> a);
double MurphyElasticCollisionCosX(double eps, std::vector<double> a);


// For Diffusion
double CalculateMeanRiRj(Particles *particle_list, int axis1, int axis2, int type);
double CalculateMeanRiVj(Particles *particle_list, int axis1, int axis2, int type);
