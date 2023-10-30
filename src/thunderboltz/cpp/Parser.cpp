void ReadInput(std::string inputFileName, InputData *SimulationParam)
{
	SimulationParam->E=0.0;
	SimulationParam->Emag=0.0;
	SimulationParam->Efreq=0.0;
	SimulationParam->B[0]=0.0;
	SimulationParam->B[1]=0.0;
	SimulationParam->B[2]=0.0;
    SimulationParam->MEM=0.0;
    SimulationParam->Slurm=0;
    SimulationParam->OS=100;
    SimulationParam->CR=0;
    SimulationParam->NSigmaMaxSamples=1000;
	std::vector<std::string> CrossSecFile,CollType;
	std::vector<int> part1, part2, prod1, prod2;
	std::vector<double> DelE;
    std::vector<std::string> velocityLoadFilenames;
	int dumpstart, dumpstride, dumptype;
    int extrapolate = 0;

	std::fstream file(inputFileName);

	while(!file.eof())
	{
		int p1, p2, pr1, pr2;
		std::string dcsf, csf,ctyp,Index;
		file >> Index; 
        DifferentialCrossSection dcs;
        dcs.modelName = "Default"; // Isotropic for I(eps,chi), one takes-all for EESD

		if(Index=="E"){
			double E;
			file >> E; 
			printf("reading electric field: %e (V/m)\n", E);
			SimulationParam->E=E;
			SimulationParam->Emag=E;
		}
		if(Index=="ET"){
			double f;
			file >> f; 
			printf("reading electric frequency: %e (Hz)\n", f);
			SimulationParam->Efreq=f;
		}
		if(Index=="B"){
			double Bx,By,Bz;
			file >> Bx >> By >> Bz ; 
			printf("reading magnetic field: %e %e %e (Tesla)\n", Bx, By, Bz);
			SimulationParam->B[0]=Bx;
			SimulationParam->B[1]=By;
			SimulationParam->B[2]=Bz;
		}
		if(Index=="L"){
			double L;
			file >> L; 
			printf("reading box length: %e (m)\n", L);
			SimulationParam->BOX_LENGTH=L;
		}

		if(Index=="VS"){
			int vsmsamples;
			file >> vsmsamples; 
			printf("reading vsigma max samples: %d\n", vsmsamples);
			SimulationParam->NSigmaMaxSamples=vsmsamples;
		}

		if(Index=="SP"){
			int NumTypes;
			file >> NumTypes; 
			printf("reading number of particle types: %d\n", NumTypes);
			SimulationParam->NUM_TYPES=NumTypes;
            // Initialize LV vector with num type
            for (int i = 0; i < NumTypes; i++) {
                SimulationParam->velocityLoadFilenames.push_back("none");
            }
		}

		if(Index=="NP"){
			int NumTypes;
            unsigned long Num;
			NumTypes=SimulationParam->NUM_TYPES;

			for(int i=0;i<NumTypes;i++)
			{
				file >> Num;
				SimulationParam->Nparticles.push_back(Num);
				printf("reading number of particles of type %i: %lu\n", i, Num);
			}
		}
		if(Index=="MP"){
			int NumTypes;
			double Mass;
			NumTypes=SimulationParam->NUM_TYPES;
			SimulationParam->Mparticles = new double[NumTypes];

			for(int i=0;i<NumTypes;i++)
			{
				file >> Mass;
				SimulationParam->Mparticles[i]=Mass;
				printf("reading mass of particle type %i: %f (amu)\n",i, Mass); // ,SimulationParam->Mparticles[i]);
			}
		}

		if(Index=="TP"){
			int NumTypes;
			double Temp;
			NumTypes=SimulationParam->NUM_TYPES;
			SimulationParam->Tparticles= new double[NumTypes];
			for(int i=0;i<NumTypes;i++)
			{
				file >> Temp;
				SimulationParam->Tparticles[i]=Temp;
				printf("reading temperature of particle type %i: %f (eV)\n",i, Temp); // ,SimulationParam->Tparticles[i]);
			}
		}

		if(Index=="VV"){
			int NumTypes;
			double flow;
			NumTypes=SimulationParam->NUM_TYPES;
			SimulationParam->Vparticles= new double[NumTypes];
			for(int i=0;i<NumTypes;i++)
			{
				file >> flow;
				SimulationParam->Vparticles[i]=flow;
				printf("reading flow velocity of particle type %i: %f (eV)\n",i, flow); // ,SimulationParam->Vparticles[i]);
			}
		}

		if(Index=="QP"){
			int NumTypes;
			double Charge;
			NumTypes=SimulationParam->NUM_TYPES;
			SimulationParam->Qparticles= new double[NumTypes];

			for(int i=0;i<NumTypes;i++)
			{
				file >> Charge;
				SimulationParam->Qparticles[i]=Charge;
				printf("reading charge of particle type %i: %f (elementary charge)\n", i, Charge);
			}
		}

		if(Index=="NS"){
			int NS;
			file >> NS; 
			printf("reading number of time steps: %d\n", NS);
			SimulationParam->NUM_TS=NS;
		}
		if(Index=="DT"){
			double stepsize;
			file >> stepsize; 
			printf("reading time step dt: %e\n", stepsize);
			SimulationParam->DT=stepsize;
		}
		if(Index=="CS"){

			double DeltaE;
			file >>csf>>p1>>p2>>ctyp>>DeltaE>>pr1>>pr2; 
			part1.push_back(p1);
			part2.push_back(p2);
			printf("found %s file for particle %i with %i in indeck", csf.c_str(), p1,p2);
			CrossSecFile.push_back(csf);
			CollType.push_back(ctyp);
			DelE.push_back(DeltaE*QE);
			prod1.push_back(pr1);
			prod2.push_back(pr2);
            // Read optional differential parameters
            std::string line;
            std::getline(file, line);
            std::stringstream ss(line);
            if (line.find_first_not_of(" ") != std::string::npos) {
                std::stringstream ss(line);
                ss >> dcs.modelName;
                printf(", differential model found: %s", dcs.modelName.c_str());
                // Number of DCS parameters
                // int np;
                // ss >> np;
                std::string tok;
                while (ss >> tok) {
                    dcs.params.push_back(std::stof(tok));
                }
                // for (int i=0; i<np; i++) {
                    // ss >> tok;
                    // dcs.params.push_back(std::stof(tok));
                // }
            }
            printf("\n");
            SimulationParam->DCS.push_back(dcs);
		}

		if(Index=="FV")
		{
			file >> dumpstart >> dumpstride >> dumptype;
		    std::vector<int> dumpInfo = {dumpstart, dumpstride, dumptype};
            SimulationParam->DumpVDFParams.push_back(dumpInfo);

			printf("reading vdf dump spec: %d %d %d\n", dumpstart, dumpstride, dumptype);
		}

        if(Index=="OS") {
            // Add option for output stride
            file >> SimulationParam->OS;
            printf("reading output stride spec: %d\n", SimulationParam->OS);
        }

        if(Index=="SE") {
            // Add option for auto dump at the end of a SLURM run
            file >> SimulationParam->Slurm;
            printf("reading auto dump spec: %d\n", SimulationParam->Slurm);
        }

        if(Index=="CR") {
            file >> SimulationParam->CR;
            printf("reading carry remainder spec: %d\n", SimulationParam->CR);
        }

        if(Index=="EX")
        {
            file >> extrapolate;
            printf("reading extrapolation option: %d\n", extrapolate);
        }

        if(Index=="MEM") {
            file >> SimulationParam->MEM;
            printf("reading memory allocation: %e (GB)\n", SimulationParam->MEM);
        }

        if(Index=="LV") {
            std::string fname;
            int fid;
            file >> fname >> fid;
            SimulationParam->velocityLoadFilenames[fid] = fname;
            printf("reading velocity init data %s for particle type %d\n",
                fname.c_str(), fid);
            // Get the line count of this file and update NP for this species
            unsigned long NP = 0;
            std::fstream f(fname);
            std::string s;
            while (std::getline(f, s)) ++NP;

            f.close();

            SimulationParam->Nparticles[fid] = NP;
            printf("Updating particle type %d count to %lu\n", fid, NP);
        }
	}
	
	SimulationParam->collisionFilename=CrossSecFile;
	SimulationParam->particle1=part1;
	SimulationParam->particle2=part2;
	SimulationParam->product1=prod1;
	SimulationParam->product2=prod2;
	SimulationParam->collisionType=CollType;
	SimulationParam->DeltaE=DelE;
    SimulationParam->Extrapolate=extrapolate;
}
