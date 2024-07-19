/*0D DSMC Code*/
#include <unistd.h>
#include "DSMC0D.h"
#include "Collision.cpp"
#include "Parser.cpp"
#include <typeinfo>

/* --------- main -------------*/
int main(int argc, char *argv[])
{
    printf("RandMax %i\n", int(RAND_MAX));
    /*variables to hold species data ions*/
    Particles AllParticle;

    /*Test: Read Input*/
    InputData Parameters;

    if (argc > 1)
    {
        std::string arg1(argv[1]);
        ReadInput(arg1,&Parameters);
    }
    else
    {
        ReadInput("N2.in",&Parameters);
    }
    CrossSections AllCrossSections;
    /*Read Cross Sections From Input*/
    ReadCrossSection(&Parameters , &AllCrossSections);
    /*load particles*/
    InitializeParticles(&AllParticle,&Parameters);

    //Calculate square of B field
    double Bsq=Parameters.B[0]*Parameters.B[0]+Parameters.B[1]*Parameters.B[1]+Parameters.B[2]*Parameters.B[2];

    printf("------------------------------------\n");

    // Solve for parameters and print necessary output
    ReportStatus(&AllParticle, &AllCrossSections, &Parameters, 0);

    for (int i=1; i <= Parameters.NUM_TS; i++)
    {

        if (Parameters.Efreq!=0.0)
        {
            UpdateEField(&Parameters,i);
        }

        // Update displacement analytically, only track the
        // first particle type for now.
        UpdateDisplacement(&AllParticle,&Parameters,0);

        //Select Boris push if magnetic field is present.
        if(Bsq!=0.0) {
            AccelerateBoris(&AllParticle,&Parameters);
        } else {
            Accelerate(&AllParticle,&Parameters);
        }

        // Perform collision routines
        SampleVSigmaMax(&AllParticle,&AllCrossSections,Parameters.NSigmaMaxSamples);
        DSMCcollideNTC(&AllParticle,&AllCrossSections,&Parameters,i);

        // Solve for parameters and print necessary output
        ReportStatus(&AllParticle, &AllCrossSections, &Parameters, i);
    }
    return 0;
}

/***** HELPER FUNCTIONS *********************************************************/
double rnd()
{
	static std::random_device rd;
	static std::mt19937_64 generator(rd());
	static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
	double rnd  = uniform_dist(generator);
    return rnd;
}

/* random number generator for now using built-in but this is not adequate for real simulations*/
// double rnd()
// {
//     // TODO: rand() includes 0, so this can return 0, investigate if this is an issue
//     return rand()/(double)RAND_MAX;
// }

void ReportStatus(Particles *AllParticle, CrossSections *AllCrossSections, InputData *Parameters, int step)
{
    //Variables for output
    double Ke_electrons, ME_electrons;
    double Mii,Kii,Txii,Tyii,Tzii,Vxii,Vyii,Vzii,Rxi0,Ryi0,Rzi0,Nii;
    double XX, YY, ZZ, XY, XZ, YZ;
    double XVX, XVY, XVZ, YVX, YVY, YVZ, ZVX, ZVY, ZVZ;

    //Dump Particles
    DumpParticles(AllParticle, Parameters, step);

    // Quantities to dump at banner stride, default is 100
    if(step % Parameters->OS == 0)
    {
        // Calculate displacements only for particle type 0
        Rxi0 = CalculateMeanDisplacement(AllParticle, 0, 0);
        Ryi0 = CalculateMeanDisplacement(AllParticle, 1, 0);
        Rzi0 = CalculateMeanDisplacement(AllParticle, 2, 0);

        // Calculate position correlations			
	    XX=CalculateMeanRiRj(AllParticle,0,0,0);
        YY=CalculateMeanRiRj(AllParticle,1,1,0);
        ZZ=CalculateMeanRiRj(AllParticle,2,2,0);
        XY=CalculateMeanRiRj(AllParticle,0,1,0);
        XZ=CalculateMeanRiRj(AllParticle,0,2,0);
		YZ=CalculateMeanRiRj(AllParticle,1,2,0);

        // Calculation position/velocity correlations
	    XVX=CalculateMeanRiVj(AllParticle,0,0,0);
	    XVY=CalculateMeanRiVj(AllParticle,0,1,0);
	    XVZ=CalculateMeanRiVj(AllParticle,0,2,0);
	    YVX=CalculateMeanRiVj(AllParticle,1,0,0);
	    YVY=CalculateMeanRiVj(AllParticle,1,1,0);
	    YVZ=CalculateMeanRiVj(AllParticle,1,2,0);
	    ZVX=CalculateMeanRiVj(AllParticle,2,0,0);
	    ZVY=CalculateMeanRiVj(AllParticle,2,1,0);
	    ZVZ=CalculateMeanRiVj(AllParticle,2,2,0);

        //Dump Fluid Quantities
        for(int ii=0;ii<Parameters->NUM_TYPES;ii++)
        {
            std::fstream outfile;

            outfile.open("Particle_Type_"+std::to_string(ii)+".dat", std::ios_base::app);
            Kii=ComputeKE(AllParticle,ii);
            Mii=Kii/Parameters->Nparticles[ii];
            Nii=Parameters->Nparticles[ii]/pow(Parameters->BOX_LENGTH,3);
            Vxii=ComputeVelocityMoment(AllParticle,ii,0);
            Vyii=ComputeVelocityMoment(AllParticle,ii,1);
            Vzii=ComputeVelocityMoment(AllParticle,ii,2);
            Txii = ComputeTemperature(AllParticle,Parameters,ii,0);
            Tyii = ComputeTemperature(AllParticle,Parameters,ii,1);
            Tzii = ComputeTemperature(AllParticle,Parameters,ii,2);
            if(step==0) {
                outfile << "t , Ki , Mi , Ni , Vxi , Vyi , Vzi , Txi , Tyi , Tzi";
                if (ii==0) { outfile << " , Rxi , Ryi , Rzi, XX, YY, ZZ, XY, XZ, YZ, XVX, XVY, XVZ, YVX, YVY, YVZ, ZVX, ZVY, ZVZ"; }
                outfile << "\n";
            }
            outfile << step*Parameters->DT<<" , "<<Kii <<" , "<< Mii <<" , "<< Nii <<" , "<<Vxii<<" , "<<Vyii<<" , "<<Vzii<<" , "<<Txii<<" , "<<Tyii<<" , "<<Tzii;

            // If writing to electron file
            if (ii==0) {
              outfile << " , "
              << std::setprecision(9) << Rxi0 << " , " << std::setprecision(9) << Ryi0 << " , " << std::setprecision(9) << Rzi0<< " , "
              << std::setprecision(9) << XX<< " , " << std::setprecision(9) << YY<< " , " << std::setprecision(9) << ZZ<< " , "
              << std::setprecision(9) << XY<< " , " << std::setprecision(9) << XZ<< " , " << std::setprecision(9) << YZ<< " , "
              << std::setprecision(9) << XVX<< " , " << std::setprecision(9) << XVY<< " , " << std::setprecision(9) << XVZ<< " , "
              << std::setprecision(9) << YVX<< " , " << std::setprecision(9) << YVY<< " , " << std::setprecision(9) << YVZ<< " , "
              << std::setprecision(9) << ZVX<< " , " << std::setprecision(9) << ZVY<< " , " << std::setprecision(9) << ZVZ;
            }
            outfile <<"\n";
        }

        //Banner Output
        Ke_electrons=ComputeKE(AllParticle,0);
        ME_electrons=ComputeKE(AllParticle,0)/Parameters->Nparticles[0];
        printf("%i, t, %e, KEe, %f, MEe, %.9e, E, %f, NP, ",
            step, step*Parameters->DT, Ke_electrons, ME_electrons,
            Parameters->E);
        for (int type = 0; type < Parameters->NUM_TYPES-1; type++) {
            printf("%lu ", AllParticle->typeList[type].particleIdentifier.size());
        }
        printf("%lu\n", AllParticle->typeList[Parameters->NUM_TYPES-1].particleIdentifier.size());

        fflush(stdout);

        //Dump Reaction Counts
        DumpReactionCounts(AllCrossSections,step);

        // Optionally run slurm exit sequence
        if (Parameters->Slurm) {
            // Check if near the end of run.
            if (SlurmExit(step, Parameters->NUM_TS)) {
                DumpParticles(AllParticle,Parameters,step);
            }
        }
    }
}
/* samples random velocity from Maxwellian distribution*/
double SampleVel(double T, double mass)
{
    double v_th = sqrt(2*K*T/mass);
    return v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
}

void InitializeParticles(Particles *particle_list, InputData *SimulationParam)
{
    // Determine total number of particles
    unsigned long NUM_PART=MEM_FAC*std::accumulate(SimulationParam->Nparticles.begin(),SimulationParam->Nparticles.end(),0.0);//NUM_ELECTRONS+NUM_IONS+NUM_NEUTRALS; //
    // Estimate of memory per particle
    printf("Memory per particle: %lu bytes\n", sizeof(Particle));
    // Memory for all particles in GB (base 2)
    printf("Memory for all initial particles: %.3e GB\n", NUM_PART*sizeof(Particle)/(MEM_FAC*pow(2, 30)));

    if (SimulationParam->MEM > 0) {
        // Allocate based on available memory
        double NP = SimulationParam->MEM*pow(2, 30)/sizeof(Particle);
        // Ensure no particle count overflow
        const unsigned long MAX_LONG = -1;
        assert(NP < MAX_LONG);
        NUM_PART = floor(NP);
    }

    printf("Allocating space for %lu particles\n", NUM_PART);
    fflush(stdout);
    particle_list->np_alloc=NUM_PART;
    particle_list->part = new Particle[NUM_PART];
    particle_list->typeList = new ParticleList[SimulationParam->NUM_TYPES];
    particle_list->np=0;


    for (int i=0;i<SimulationParam->NUM_TYPES; i++)
    {
        std::string velFile = SimulationParam->velocityLoadFilenames[i];
        double qp = SimulationParam->Qparticles[i]*QE;
        double mp = SimulationParam->Mparticles[i]*AMU;
        double tp = SimulationParam->Tparticles[i]*EV_TO_K;
        double fp = SimulationParam->Vparticles[i]*QE;

        if (velFile != "none") {
            // Read velocities from a file
            std::fstream file(velFile);
            double vx, vy, vz = 0;
            std::string line;
            std::getline(file, line);
            while (!file.eof()) {
                // Filter out commas 
                std::replace(line.begin(), line.end(), ',', ' ');
                std::stringstream ss(line);
                ss >> vx >> vy >> vz;
                AddParticle(SimulationParam,particle_list,0,0,0,vx,vy,vz,qp,mp,i,0);
                // Read next line
                std::getline(file, line);
            }
            file.close();

        } else {
            // Generate velocities from temperature and flow velocity settings
            double FlowShift = 0;
            //convert flow in eV to m/s
            FlowShift = sqrt(2.0*fp/(mp));

            for (unsigned long p=0;p<SimulationParam->Nparticles[i];p++)
            {
                double vx = SampleVel(tp, mp);
                double vy = SampleVel(tp, mp);
                double vz = SampleVel(tp, mp)+FlowShift;
                // Create particle at 0, 0, 0 position with sampled velocities
                AddParticle(SimulationParam,particle_list,0,0,0,vx,vy,vz,qp,mp,i,0);
            }
        }
    }
}

/*adds new particle to the species, returns pointer to the newly added data*/
void AddParticle(InputData *SimulationParam, Particles *particle_list,
                 double rx, double ry, double rz, double vx, double vy, double vz,
                 double charge, double mass, int type, int step)
{
    /*abort the simulation if we ran out of space to store this particle*/
    if (particle_list->np > particle_list->np_alloc-1)
    {
        DumpParticles(particle_list,SimulationParam,step);
        printf("Too many particles!\n"); exit(-1);
    }

    /*store position and velocity of this particle*/
    particle_list->part[particle_list->np].r[0] = rx;
    particle_list->part[particle_list->np].r[1] = ry;
    particle_list->part[particle_list->np].r[2] = rz;
    particle_list->part[particle_list->np].v[0] = vx;
    particle_list->part[particle_list->np].v[1] = vy;
    particle_list->part[particle_list->np].v[2] = vz;
    particle_list->part[particle_list->np].mass = mass;
    particle_list->part[particle_list->np].charge = charge;
    particle_list->part[particle_list->np].type = type;

    //add to type list
    particle_list->typeList[type].particleIdentifier.push_back(particle_list->np);
    /*increment particle counter*/
    particle_list->np++;
}

/* Accelerate Species*/
void Accelerate(Particles *particle_list, InputData *SimulationParam)
{
    // Declare q/m
    double qm = 0;
    // Loop over particles
    for (unsigned long p=0;p<particle_list->np;p++)
    {
        // Grab pointer to this particle
        Particle *part = &particle_list->part[p];

        qm=part->charge / part->mass;
        // Advance velocity
        part->v[2] += SimulationParam->DT*qm*SimulationParam->E;
    }
}

void AccelerateBoris(Particles *particle_list, InputData *SimulationParam)
{
    /*precompute q/m*/
    double qm = 0;
    double t[3], s[3], v_minus[3], v_minus_cross_t[3], v_prime[3], v_prime_cross_s[3], v_plus[3];
    double tsq;
    /*loop over particles*/
    for (unsigned long p=0;p<particle_list->np;p++)
    {
        /*grab pointer to this particle*/
        Particle *part = &particle_list->part[p];

        qm=part->charge / part->mass;
        /*advance velocity*/

        /*t vector*/
        for (int dim=0;dim<3;dim++)
        {
            t[dim] = qm*SimulationParam->B[dim]*0.5*SimulationParam->DT;
        }

        tsq=t[0]*t[0]+t[1]*t[1]+t[2]*t[2];

        /*s vector*/
        for (int dim=0;dim<3;dim++)
        {
            s[dim] = 2*t[dim]/(1+tsq);
        }

        /*v minus*/
        v_minus[0] = part->v[0];
        v_minus[1] = part->v[1];
        v_minus[2] = part->v[2] + qm*SimulationParam->E*0.5*SimulationParam->DT;

        /*v_minus cross t*/
        v_minus_cross_t[0]=v_minus[1]*t[2]-v_minus[2]*t[1];
        v_minus_cross_t[1]=-v_minus[0]*t[2]+v_minus[2]*t[0];
        v_minus_cross_t[2]=v_minus[0]*t[1]-v_minus[1]*t[0];

        /*v prime*/
        for (int dim=0;dim<3;dim++)
        {
            v_prime[dim] = v_minus[dim] + v_minus_cross_t[dim];
        }

        /*v_prime cross s*/
        v_prime_cross_s[0]=v_prime[1]*s[2]-v_prime[2]*s[1];
        v_prime_cross_s[1]=-v_prime[0]*s[2]+v_prime[2]*s[0];
        v_prime_cross_s[2]=v_prime[0]*s[1]-v_prime[1]*s[0];

        //========

        for (int dim=0;dim<3;dim++)
        {
            v_plus[dim] = v_minus[dim] + v_prime_cross_s[dim];
        }

        part->v[0] = v_plus[0];
        part->v[1] = v_plus[1];
        part->v[2] = v_plus[2] + qm*SimulationParam->E*0.5*SimulationParam->DT;
    }
}

void UpdateDisplacement(Particles *particles,InputData *SimulationParam, int type) {
    double dt = SimulationParam->DT;
    double E = SimulationParam->E;
    for (unsigned long p=0;p<particles->np;p++) {
        if (particles->part[p].type == type) {
            Particle *pt = &particles->part[p];
            // Update position based on velocity and acceleration
            pt->r[0] += pt->v[0]*dt;
            pt->r[1] += pt->v[1]*dt;
            pt->r[2] += pt->v[2]*dt + 0.5*pt->charge*E/pt->mass*dt*dt;
        }
    }
}

double CalculateMeanDisplacement(Particles *particle_list, int axis, int type) {
    double r,N;
    N=0;
    r=0;
    for (unsigned long p=0;p<particle_list->np;p++) {
        if (particle_list->part[p].type == type) {
            r += particle_list->part[p].r[axis];
            N += 1.0;
        }
    }
    r=r/N;
    return r;
}

// For bulk diffusion
double CalculateMeanRiRj(Particles *particle_list, int axis1, int axis2, int type)
{
    double xi,xj,xy,N;
    N=0.0;
    xy=0.0;
	for (unsigned long p=0;p<particle_list->np;p++)
	{
        if (particle_list->part[p].type == type)
        {
			xi=particle_list->part[p].r[axis1];
			xj=particle_list->part[p].r[axis2];
			xy+=xi*xj;
			N+=1.0;
        }
	}
	xy=xy/(N);
	return xy;
}

// For flux diffusion
double CalculateMeanRiVj(Particles *particle_list, int axis1, int axis2, int type)
{
    double xi, vj, xvy, N;
    N=0.0;
    xvy=0.0;
	for (unsigned long p=0;p<particle_list->np;p++)
	{
        if (particle_list->part[p].type == type)
        {
			xi=particle_list->part[p].r[axis1];
			vj=particle_list->part[p].v[axis2];
			xvy+=xi*vj;
			N+=1.0;
        }
	}
	xvy=xvy/N;
	return xvy;
}

/*computes species kinetic energy in electron volt*/
double ComputeKE(Particles *particle_list, int type)
{
    double ke=0;
    double vmagsq=0;
    for (unsigned long p=0;p<particle_list->np;p++)
        if (particle_list->part[p].type == type)
        {
            vmagsq = (particle_list->part[p].v[0]*particle_list->part[p].v[0]
                    + particle_list->part[p].v[1]*particle_list->part[p].v[1]
                    + particle_list->part[p].v[2]*particle_list->part[p].v[2]);
            /*we now have sum of v^2, multiply by 0.5*mass*/
            ke += 0.5*vmagsq*particle_list->part[p].mass;
        }
    /*convert to electron volts, 1eV=QE joules*/
    ke /= QE;
    return ke;
}


double ComputeMeanE(Particles *particle_list, InputData *SimulationParam, int type)
{
    double me=0;
    double vmagsq=0;
    double vsq=0;
    for (unsigned long p=0;p<particle_list->np;p++)
        if (particle_list->part[p].type == type)
        {
            vmagsq = (particle_list->part[p].v[0]*particle_list->part[p].v[0]
                    + particle_list->part[p].v[1]*particle_list->part[p].v[1]
                    + particle_list->part[p].v[2]*particle_list->part[p].v[2]);
            vsq += vmagsq;
        }
    /*we now have sum of v^2, multiply by 0.5*mass*/
    /*convert to electron volts, 1eV=QE joules*/
    me = vsq/SimulationParam->Nparticles[type]*SimulationParam->Mparticles[type]*AMU/(2*SimulationParam->Qparticles[type]*QE);
    return me/QE;
}

/*computes species average velocity in one coordinate axis*/
double ComputeVelocityMoment(Particles *particle_list, int type, int axis)
{
    double Vmoment=0;
    double n=0;
    for (unsigned long p=0;p<particle_list->np;p++)
        if (particle_list->part[p].type == type)
        {
            Vmoment+=particle_list->part[p].v[axis];
            n+=1;
        }
    Vmoment=Vmoment/n;
    return Vmoment;
}

/*Reads Cross Section Data*/
void ReadCrossSection(InputData *SimulationParam, CrossSections *CrossSecList)
{

    int numcrosssec;
    /* number of cross sections is the same as the size of the lists particle 1 or particle 2*/
    numcrosssec=SimulationParam->particle1.size();
    //
    /* Create Arrays for Cross Sections and VSigmaMax. Initialize Vsigmamax later after particles vdfs are initialized*/
    CrossSecList->csec = new CrossSec[numcrosssec];
    CrossSecList->VSigmaMax = new double[numcrosssec];
    // Create array for Npair remainders and initialize to 0.
    CrossSecList->remainders = new double[numcrosssec];
    for (int i = 0; i < numcrosssec; i++) CrossSecList->remainders[i] = 0;

    //Read in collision data file names and collision particle types from SimulationParameters
    std::vector<std::string> filename,collType;
    std::vector<int> type1, type2, product1,product2;
    std::vector<double> DeltaE;

    DeltaE=SimulationParam->DeltaE;
    filename=SimulationParam->collisionFilename;
    type1=SimulationParam->particle1;
    type2=SimulationParam->particle2;
    product1=SimulationParam->product1;
    product2=SimulationParam->product2;
    collType=SimulationParam->collisionType;

    std::vector<double> sig;            /*array holding cross sections, energies*/
    std::vector<double> ener;
    CrossSecList->NCrossSec=filename.size();

    for (unsigned int i = 0; i < filename.size(); i++) 
    {
        std::fstream file(filename[i]);
        printf("reading %s\n", filename[i].c_str());

        //Make sure energy and cross sec vectors are clear
        sig.clear();
        ener.clear();

        while(!file.eof())
        {
            double a, b;
            file >> a >> b; // extracts 2 floating point values seperated by whitespace
            // printf("%d %.13e %.13e\n", i, a, b);
            ener.push_back(a);
            sig.push_back(b);
        }

        CrossSecList->csec[i].sigma =sig;
        CrossSecList->csec[i].energy=ener;
        CrossSecList->csec[i].ParticleType1=type1[i];
        CrossSecList->csec[i].ParticleType2=type2[i];
        CrossSecList->csec[i].ProductType1=product1[i];
        CrossSecList->csec[i].ProductType2=product2[i];

        CrossSecList->csec[i].collisionType=collType[i];
        CrossSecList->csec[i].DeltaE=DeltaE[i];
        CrossSecList->csec[i].count=0;
    }
}

void DSMCcollideNTC(Particles *particle_list, CrossSections *CrossSecList, InputData *SimulationParam, int step)
{
    int numcollisiontype;
    int particleType1, particleType2,product2Type,product1Type;
    numcollisiontype=CrossSecList->NCrossSec;

    for(int i=0;i<numcollisiontype; i++)
    {
        if(quiet==1) printf("==========================\n");
        if(quiet==1) printf("Collision Diagnostics\n");
        if(quiet==1) printf("==========================\n");
        if(quiet==1) printf("collision type %i \n",i+1);

        std::string collisionType = CrossSecList->csec[i].collisionType;
        // printf("%s\n", collisionType.c_str());
        particleType1=CrossSecList->csec[i].ParticleType1;
        particleType2=CrossSecList->csec[i].ParticleType2;
        product2Type = CrossSecList->csec[i].ProductType2;
        product1Type = CrossSecList->csec[i].ProductType1;
        if(quiet==1) printf("colliding particle type %i with type %i \n",particleType1,particleType2);

        // Store type list aliases
        std::vector<unsigned long> &partType1IDs = particle_list->typeList[particleType1].particleIdentifier;
        std::vector<unsigned long> &partType2IDs = particle_list->typeList[particleType2].particleIdentifier;
        std::vector<unsigned long> &prodType1IDs = particle_list->typeList[product1Type].particleIdentifier;
        std::vector<unsigned long> &prodType2IDs = particle_list->typeList[product2Type].particleIdentifier;

        int N1=partType1IDs.size();
        int N2=partType2IDs.size();
        if(quiet==1) printf("N1 = %i, N2 = %i \n", N1, N2);

        double delta12 = 1.0; if(particleType1==particleType2) delta12 = 0.5;

        // Calculate Npairs as a float and save the remainder before rounding
        double NP = double(N1)*double(N2) * delta12/pow(SimulationParam->BOX_LENGTH,3.0)
                    * CrossSecList->VSigmaMax[i] * SimulationParam->DT;
        // Add previous remainder to current Npairs (per Bird NTC to maintain correct sampling distribution)
        double &rem = CrossSecList->remainders[i];
        double NP_rem = NP + rem;
        unsigned long Npairs = floor(NP_rem);

        // Update the new remainder if CR param is set
        if (SimulationParam->CR) rem = NP_rem - Npairs;

        if(quiet==1) printf("Npairs %lu \n", Npairs);

        for (unsigned long j=0; j<Npairs; j++) {
            //Select two random particles for the test collision pair
            unsigned long ParticleID1; // Store the ID so that we don't need to search for it later
            unsigned long RandomParticleType1 = particle_list->RandOfType(particleType1, &ParticleID1);
            unsigned long ParticleID2;
            unsigned long RandomParticleType2 = particle_list->RandOfType(particleType2, &ParticleID2);

            // Compute relative energies in eV
            double E12 = particle_list->Erel(RandomParticleType1,RandomParticleType2);
            // Compute relative velocity in m/s
            double V12 = particle_list->Vrel(RandomParticleType1,RandomParticleType2);
            // Fetch total cross section in m^2
            double S12 = CrossSecList->csec[i].SigmaOfE(E12,SimulationParam->Extrapolate);

            // Probability that this pair will undergo this process `collisionType`
            double Ppair=V12*S12/CrossSecList->VSigmaMax[i];

            if(quiet==1) printf("probability P%lu = %f \n",j,Ppair);

            if (Ppair>rnd()) 
            {
                if(quiet==1) printf("collision pass\n");

                if (collisionType=="Elastic") {
                    // Run elastic scattering subroutine, update both particle velocities
                    ElasticCollision(particle_list,RandomParticleType1,RandomParticleType2,false,SimulationParam->DCS[i]);

                } else if (collisionType=="ElasticFixedParticle2") {
                    // Run elastic scattering subroutine, update both particle velocities
                    ElasticCollision(particle_list,RandomParticleType1,RandomParticleType2,true,SimulationParam->DCS[i]);

                } else if(collisionType=="Inelastic") {
                    InelasticCollision(particle_list,RandomParticleType1,RandomParticleType2, CrossSecList->csec[i].DeltaE,false);

                } else if (collisionType=="InelasticChangeParticle2") {
                    InelasticCollision(particle_list,RandomParticleType1,RandomParticleType2 ,CrossSecList->csec[i].DeltaE, false);

                    //Change the type of particle 2 to the product
                    particle_list->part[RandomParticleType2].type=product2Type;
                    //erase element in one type list
                    partType2IDs.erase(partType2IDs.begin()+ParticleID2);
                    //add element in the other type list
                    prodType2IDs.push_back(RandomParticleType2);

                    //Update N particles
                    SimulationParam->Nparticles[particleType2]=partType2IDs.size();
                    SimulationParam->Nparticles[product2Type]=prodType2IDs.size();

                } else if(collisionType=="InelasticChangeParticles") {
                    InelasticCollision(particle_list,RandomParticleType1,RandomParticleType2 ,CrossSecList->csec[i].DeltaE, false);
                    //Change the type of particle 2 to the product
                    particle_list->part[RandomParticleType2].type=product2Type;
                    //erase element in one type list
                    partType2IDs.erase(partType2IDs.begin()+ParticleID2);
                    //add element in the other type list
                    prodType2IDs.push_back(RandomParticleType2);

                    //Update N particles
                    SimulationParam->Nparticles[particleType2]=partType2IDs.size();
                    SimulationParam->Nparticles[product2Type]=prodType2IDs.size();
                    //Change the type of particle 1 to the product
                    particle_list->part[RandomParticleType1].type=product1Type;
                    //erase element in one type list
                    partType1IDs.erase(partType1IDs.begin()+ParticleID1);
                    //add element in the other type list
                    prodType1IDs.push_back(RandomParticleType1);

                    //Update N particles
                    SimulationParam->Nparticles[particleType1]=partType1IDs.size();
                    SimulationParam->Nparticles[product1Type]=prodType1IDs.size();

                } else if (collisionType=="N2Dissociation") {
                    //Assume reaction is in order e + N2 -> e + N + N

                    //step 1: inelastic collision with electron
                    //-----------------------------------------
                    InelasticCollision(particle_list,RandomParticleType1,RandomParticleType2, CrossSecList->csec[i].DeltaE, false);

                    //step 2: remove product N2^* 
                    //-----------------------------------------
                    // Store position of scattered electron.
                    double XReac0=particle_list->part[RandomParticleType2].r[0];
                    double XReac1=particle_list->part[RandomParticleType2].r[1];
                    double XReac2=particle_list->part[RandomParticleType2].r[2];
                    //Save KE of product before removing
                    double VProd0=particle_list->part[RandomParticleType2].v[0];
                    double VProd1=particle_list->part[RandomParticleType2].v[1];
                    double VProd2=particle_list->part[RandomParticleType2].v[2];
                    double MProd2=particle_list->part[RandomParticleType2].mass;

                    //erase element in one type list
                    partType2IDs.erase(partType2IDs.begin()+ParticleID2);
                    //Update N particles
                    SimulationParam->Nparticles[particleType2]=partType2IDs.size();

                    //step 3: create N + N with KE of N2^* + energy of dissociation
                    //-----------------------------------------
                    //Calculate random direction for splitting
                    double eps = rnd()* 2*PI;
                    double cosX = 2.0*rnd() - 1.0;
                    double sinX = sqrt(1.0 - cosX*cosX);

                    //Calculate velocity kick due to dissociation
                    double DeltaEDiss=0.9*QE; // 0.9eV in J
                    double VDissMag=sqrt(2.0*DeltaEDiss/MProd2);

                    //using angle and VDissMag calculate complimentary velocity kicks
                    double dvx = VDissMag*cosX;
                    double dvy = VDissMag*sinX*cos(eps);
                    double dvz = VDissMag*sinX*sin(eps);

                    //Create two new particles, for now no tracking the neutrals
                    AddParticle(SimulationParam,particle_list,XReac0,XReac1,XReac2,VProd0+dvx,VProd1+dvy,VProd2+dvz,0.0,14.0*AMU,product2Type,step);
                    AddParticle(SimulationParam,particle_list,XReac0,XReac1,XReac2,VProd0-dvx,VProd1-dvy,VProd2-dvz,0.0,14.0*AMU,product2Type,step);

                    //Update N Particles
                    SimulationParam->Nparticles[product2Type]=prodType2IDs.size();

                } else if (collisionType.find("Ionization") != std::string::npos) {
                    //Assume reaction is in order e + n -> e + n^+ + e, change particle 2 like in inelastic collision but then add electron.
                    // Save velocities before collision
                    double (&v1)[3] = particle_list->part[RandomParticleType1].v; // incident electron velocity
                    double (&v2)[3] = particle_list->part[RandomParticleType2].v; // target macroparticle velocity
                    // Save position of incident particle
                    double (&r1)[3] = particle_list->part[RandomParticleType1].r; // Incident electron position
                    // Use the exact electron mass provided
                    double m_e = particle_list->part[RandomParticleType1].mass;

                    // Relative velocities
                    double g[3];
                    vadd_scale(g, v1, v2, -1); // Subtract the two vectors
                    double gsqd = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
                    double gmag = sqrt(gsqd);

                    // Obtain lab frame euler angle ratios, avoid divide by 0.
                    double sin_theta;
                    double cos_theta;
                    double phi;
                    if (g[0] == 0) cos_theta = 0; else cos_theta = g[0]/gmag;
                    sin_theta = sqrt(1 - fmin(1, cos_theta*cos_theta));
                    if (g[1] == 0) {
                        if (g[2] > 0) phi = 0.5*PI; else phi = -0.5*PI;
                    } else {phi = atan2(g[2], g[1]);}
                    double sin_phi = sin(phi);
                    double cos_phi = cos(phi);

                    // The energy of the incident, scattered and ejected electrons (in Joules)
                    // Calculate assuming rest frame of target, and assuming mu ~= m_e
                    double ep = 0.5*m_e*gsqd;
                    double ep_s;
                    double ep_ej;
                    // Ionization threshold in Joules
                    double ep_ion = CrossSecList->csec[i].DeltaE;
                    if (SimulationParam->DCS[i].modelName == "Default") {
                        ep_ej = 0;
                    } else if (SimulationParam->DCS[i].modelName == "Equal") {
                        // Equal energy sharing, keep everythin in Joules
                        ep_ej = (ep-ep_ion)/2;
                    } else if (SimulationParam->DCS[i].modelName == "Uniform") {
                        ep_ej = rnd()*(ep-ep_ion)/2;
                    } else {
                        throw std::invalid_argument(
                            "Received unimplemented ionization energy sharing model");
                    }
                    // The energy of the primary / scattered electron
                    ep_s = ep - ep_ion - ep_ej; // These are in COM frame

                    // Assume that these energies in the COM frame are
                    // approximately the same in the target frame due to the
                    // large mass of the macroparticle

                    // Sample post collision angles for scattered electron
                    double cos_chi_s = 2*rnd() - 1;
                    double sin_chi_s = sqrt(1 - pow(cos_chi_s, 2));
                    double eta_s = 2*PI*rnd();
                    double cos_eta_s = cos(eta_s);
                    double sin_eta_s = sin(eta_s);
                    // Sample post collision angles for ejected electron
                    double cos_chi_ej = 2*rnd() - 1;
                    double sin_chi_ej = sqrt(1 - pow(cos_chi_ej, 2));
                    double eta_ej = 2*PI*rnd();
                    double cos_eta_ej = cos(eta_ej);
                    double sin_eta_ej = sin(eta_ej);
                    // Get velocity vector magnitude of g' from energy of each particle
                    double gmag_s = sqrt(2*ep_s/m_e);
                    double gmag_ej = sqrt(2*ep_ej/m_e);
                    // Build post collision unit vector g^'_T in rotated frame for ej and scat. e-
                    double gp_s[3];
                    double gp_ej[3];
                    // Convert angles to cartesian coordinates, (still a unit vector)
                    Project(gp_s, cos_chi_s, sin_chi_s, cos_eta_s, sin_eta_s);
                    Project(gp_ej, cos_chi_ej, sin_chi_ej, cos_eta_ej, sin_eta_ej);
                    if (false) {  // This stage is only necessary if the angular distribution is anisotropic
                        // Rotate velocity back to the unrotated COM frame
                        RotateFrame(gp_s, cos_theta, sin_theta, cos_phi, sin_phi);
                        RotateFrame(gp_ej, cos_theta, sin_theta, cos_phi, sin_phi);
                    }
                    // Update scattered electron velocities in lab frame and scale unit vector with magnitude gmag_s
                    vadd_scale(v1, v2, gp_s, gmag_s);

                    // Change the type of particle 2 to the product
                    particle_list->part[RandomParticleType2].type=product2Type;
                    // erase element in one type list
                    partType2IDs.erase(partType2IDs.begin()+ParticleID2);
                    // add element in the other type list
                    prodType2IDs.push_back(RandomParticleType2);

                    if (collisionType.find("NoEGen") == std::string::npos) {
                        // Create velocity for ejected electron
                        double v_ej[3];
                        // Update ejected electron velocities in lab frame and scale unit vector with magnitude gmag_s
                        vadd_scale(v_ej, v2, gp_ej, gmag_ej);
                        // Add particle with ejected electron velocities
                        AddParticle(SimulationParam, particle_list, r1[0], r1[1], r1[2], v_ej[0], v_ej[1], v_ej[2], SimulationParam->Qparticles[0]*QE,SimulationParam->Mparticles[0]*AMU,0,step);
                        // Update N particles
                        SimulationParam->Nparticles[particleType2]=partType2IDs.size();
                        SimulationParam->Nparticles[product2Type]=prodType2IDs.size();
                        // update electron size
                        SimulationParam->Nparticles[0]=particle_list->typeList[0].particleIdentifier.size();
                        // Set charge of ion
                        particle_list->part[RandomParticleType2].charge=SimulationParam->Qparticles[product2Type]*QE;
                    }
                } else if (collisionType=="InelasticFixedParticle2") {
                    InelasticCollision(particle_list,RandomParticleType1,RandomParticleType2, CrossSecList->csec[i].DeltaE, true);

                } else if (collisionType=="CX") {
                    CXCollision(particle_list,RandomParticleType1,RandomParticleType2);
                }

                CrossSecList->csec[i].count+=1; //increment interaction count
            } else {
                if(quiet==1) printf("collision fail\n");
            }
        }
    }
}

double Particles::Vrel(unsigned long partid1, unsigned long partid2)
{
    double vrx=part[partid1].v[0]-part[partid2].v[0];
    double vry=part[partid1].v[1]-part[partid2].v[1];
    double vrz=part[partid1].v[2]-part[partid2].v[2];
    double vr=pow(vrx*vrx+vry*vry+vrz*vrz,0.5);
    return vr;
}

double Particles::Erel(unsigned long partid1, unsigned long partid2)
{
    double vrx=part[partid1].v[0]-part[partid2].v[0];
    double vry=part[partid1].v[1]-part[partid2].v[1];
    double vrz=part[partid1].v[2]-part[partid2].v[2];
    double vrsq=vrx*vrx+vry*vry+vrz*vrz;
    double m1 = part[partid1].mass;
    double m2 = part[partid2].mass;
    double mu = m1*m2/(m1+m2);
    return 0.5*mu*vrsq/QE; // In eV
}

int Particles::RandOfType(int type, unsigned long *pid)
{
    // Return the 
    int rand, randpart;

    rand = int(rnd() * typeList[type].particleIdentifier.size());
    randpart = typeList[type].particleIdentifier[rand]; /*??*/
    *pid = rand; // save the index in the type list

    /* For Debugging:*/
    int testtype;
    testtype = part[randpart].type;
    if(testtype!=type) printf("mismatch: testtype = %i type = %i\n",testtype,type);

    return randpart;
}

double CrossSec::SigmaOfE(double PartEnergy, bool extrapolate)
{
    //adapted from http://www.cplusplus.com/forum/general/216928/
   int size = energy.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( PartEnergy >= energy[size - 2] )                                                 // special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while ( PartEnergy > energy[i+1] ) i++;
   }
   double xL = energy[i], yL = sigma[i], xR = energy[i+1], yR = sigma[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( PartEnergy < xL ) yR = yL;
      if ( PartEnergy > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( PartEnergy - xL );                                              // linear interpolation
}

void SampleVSigmaMax(Particles *particle_list, CrossSections *CrossSecList, int samples)
{
    /*currently assumes one cross section per species pair*/

    for(int i=0; i < CrossSecList->NCrossSec; i++)
    {
        /*identify particle types takeing part in collision*/
        int p1=CrossSecList->csec[i].ParticleType1;
        int p2=CrossSecList->csec[i].ParticleType2;

        //check to see that the particle densities are not zero
        int Np1,Np2;
        Np1=particle_list->typeList[p1].particleIdentifier.size();
        Np2=particle_list->typeList[p2].particleIdentifier.size();

        double vSigmaMax=0.0;

        //sample only if Np1 and Np2 are nonzero, otherwise vSigmaMax=0
        if(Np1*Np2!=0)
        {
        /*sample vsigma */
            for(int j = 1; j<samples; j++)
            {
                unsigned long dummy; // We do not need the type list index
                int Rp1 = particle_list->RandOfType(p1, &dummy);
                int Rp2 = particle_list->RandOfType(p2, &dummy);
                double E12 = particle_list->Erel(Rp1,Rp2); // eV
                double V12 = particle_list->Vrel(Rp1,Rp2);
                double S12 = CrossSecList->csec[i].SigmaOfE(E12,false);
                // printf("Rand particle 1: %i , Rand particle 2: %i, Vrel 1 %f, Erel %f, sigma %e\n", Rp1, Rp2,V12,E12,S12);
                if (V12*S12>vSigmaMax) vSigmaMax=V12*S12;
            }
        }
        /*save new sample vSigmaMax*/
        CrossSecList->VSigmaMax[i]=vSigmaMax;
    }
}

double ComputeTemperature(Particles *particle_list, InputData *SimulationParam,int type, int axis)
{

    double Vmoment=0;
    double n=0;
    for (unsigned long p=0;p<particle_list->np;p++)
        if (particle_list->part[p].type == type)
        {
            Vmoment+=particle_list->part[p].v[axis];
            n+=1;
        }
    Vmoment=Vmoment/n;
    
    /*compute kinetic energy of the thermal component*/
    double ke=0;
    for (unsigned long p=0;p<particle_list->np;p++)
    {
        if (particle_list->part[p].type == type)
        {
        double dv = particle_list->part[p].v[axis]-Vmoment;
        ke += dv*dv;
        }
    }

    /*we now have sum of v^2, get average kinetic energy*/
    double ke_ave;
    ke_ave=0.5*SimulationParam->Mparticles[type]*AMU*ke/SimulationParam->Nparticles[type];

    /*compute temperature in eV*/
    double T = 2*ke_ave/QE;

    return T;
}

void DumpParticles(Particles *particle_list, InputData *Parameters, int step)
{
    std::fstream outfile;
    std::vector<std::vector<int>> &vdfParams = Parameters->DumpVDFParams;
    int start, stride, type, isRequested;
    // Loop through the particle types
    for (int ptype=0; ptype < Parameters->NUM_TYPES; ptype++) {
        isRequested = 0;
        // And each FV request
        for (auto &vdfParam : vdfParams) {
            start = vdfParam[0]; stride = vdfParam[1]; type = vdfParam[2];
            if (ptype != type) continue;
            isRequested = isRequested | (step >= start && (step-start) % stride == 0);
        }

        if (isRequested) {

            outfile.open("Particle_"+std::to_string(ptype)+"_Dump_"+std::to_string(step)+".dat", std::ios_base::app);
            for (unsigned long p=0;p<particle_list->np;p++) {
                if (particle_list->part[p].type == ptype)
                {
                   outfile << particle_list->part[p].v[0] << " , "
                           << particle_list->part[p].v[1] << " , "
                           << particle_list->part[p].v[2] << "\n";
                }
            }
            outfile.close();
        }
    }
}

int SlurmExit(int i, int NS) {
    /* Return true if the run will end roughly before the next banner output,
     * or currently at last time step */
    // Store time left, string and char dumps
    unsigned int jobnum;
    char cmd[100];
    std::string g;

    // Read in jobnum
    std::fstream j("jobnum");
    std::getline(j, g);
    j >> g >> g >> g >> jobnum;
    j.close();

    // Get job time data
    std::snprintf(cmd, 100, "squeue --job %d -o %%L > job_status.txt", jobnum);
    std::system(cmd);
    std::snprintf(cmd, 100, "squeue --job %d -o %%l >> job_status.txt", jobnum);
    std::system(cmd);

    // Create file reader for job time status
    std::fstream f("job_status.txt");
    std::getline(f, g);
    std::getline(f, g);
    // Compute time left in seconds
    unsigned int time_left = ReadTime(g);
    // Read in the time so far
    std::getline(f, g);
    std::getline(f, g);
    unsigned int time_sf = ReadTime(g) - time_left;
    f.close();

    // Log time
    std::ofstream timeout("time_log.dat", std::ios_base::app); // Append mode
    timeout << i << " " << time_sf << "\n";
    timeout.close();
    // Read each timestep
    std::ifstream timein("time_log.dat");
    std::vector<int> tics;
    std::vector<int> time;
    int tint, time_int;
    while (timein >> tint >> time_int) {
        tics.push_back(tint);
        time.push_back(time_int);
    }
    timein.close();

    // Check last 2 time steps for time change
    int s = tics.size();
    if (s >= 3) {
        // Calculate the estimated time until next banner output
        double d1 = time[s-2] - time[s-3];
        double d2 = time[s-1] - time[s-2];
        double time_till = (d2 + (d2-d1));
        // Pad the estimate with a factor of 2.
        // printf("%d %e %d %d\n", i, time_till, time_left, jobnum);
        if ((i+1 == NS) | ((2*time_till > time_left) & (time_sf > 0))) {
            return 1;
        }
    }
    return 0;
}

void DumpReactionCounts(CrossSections *CrossSecList,int step)
{
    std::fstream outfile;
    outfile.open("Counts.dat", std::ios_base::app);
    outfile<< step ;

    for(int j = 0; j<CrossSecList->NCrossSec; j++)
    {
        outfile<< "," << CrossSecList->csec[j].count; 
    };
    outfile <<"\n";

}

void UpdateEField(InputData *SimulationParam, int step)
{
    double Enew;
    Enew=SimulationParam->Emag*sin(SimulationParam->Efreq*2.0*PI*step*SimulationParam->DT);
    SimulationParam->E=Enew;
}

unsigned int ReadTime(std::string line) {
    // Assume time is >= 60 seconds in the form MM:SS or HH:MM:SS, return seconds
    unsigned int H = 0;
    unsigned int M = 0;
    unsigned int S = 0;
    if (line.size() > 5) {
        sscanf(line.c_str(), "%d:%d:%d", &H, &M, &S);
    } else {
        sscanf(line.c_str(), "%d:%d", &M, &S);
    }
    return 3600*H + 60*M + S;
}
