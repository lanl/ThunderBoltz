void ElasticCollision(Particles *particle_list, unsigned long RandomParticleType1, unsigned long RandomParticleType2, bool FixedHeavyParticle, DifferentialCrossSection dcs)
{
	// Create references to colliding particle intitial velocities
    double (&v1)[3] = particle_list->part[RandomParticleType1].v;
    double (&v2)[3] = particle_list->part[RandomParticleType2].v;
	// Particle masses, reduced mass and other needed quantities.
	double m1 = particle_list->part[RandomParticleType1].mass;
	double m2 = particle_list->part[RandomParticleType2].mass;
    double inv_m1m2 = 1/(m1+m2);
    double m2_m1m2 = m2*inv_m1m2;
    double mr = m1*m2_m1m2;

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

	// Total kinetic energy of the system in COM frame
	double etrans = 0.5*gsqd*mr; // J
    // Center of mass velocity
    double w[3] = {((m1*v1[0])+(m2*v2[0])) * inv_m1m2,
                   ((m1*v1[1])+(m2*v2[1])) * inv_m1m2,
                   ((m1*v1[2])+(m2*v2[2])) * inv_m1m2};

    // Create isotropic azimuthal angle
	double eta = rnd()*2*PI;
    double sin_eta = sin(eta);
    double cos_eta = cos(eta);
    // Determine cos component of scattering angle from `modelName`
    double cos_chi = 0;
    if (dcs.modelName == "Default") { // Isotropic
        cos_chi = 2.0*rnd() - 1.0;
    } else if (dcs.modelName == "Murphy") { // Assume all of the following models are anisotropic
        cos_chi = MurphyElasticCollisionCosX(etrans/QE, dcs.params);
    } else if (dcs.modelName == "Park") {
        cos_chi = ParkElasticCollisionCosX(etrans/QE, dcs.params);
    } else {
        throw std::invalid_argument("Received unimplemented differential elastic collision model");
    }
    // Determine sin component of scattering angle
    // Sometimes cos_chi < -1.0 from precision error
    double sin_chi = sqrt(1 - fmin(1, cos_chi*cos_chi));

    // Build post collision unit vector g^'_T in rotated COM frame
    double gp[3];
    // Convert angles to cartesian coordinates
    Project(gp, cos_chi, sin_chi, cos_eta, sin_eta);
    if (dcs.modelName != "Default") {  // This stage is only necessary if the angular distribution is anisotropic
        // Rotate velocity back to the unrotated COM frame
        RotateFrame(gp, cos_theta, sin_theta, cos_phi, sin_phi);
    }
    // Set the resulting vector in the lab frame from COM velocities / reduced mass
    vadd_scale(v1, w, gp, m2_m1m2*gmag);

    // printf("END %e\n", particle_list->part[RandomParticleType1].v[2]);

    // Change second particle velocity if necessary
	if(!FixedHeavyParticle){ vadd_scale(v2, w, gp, -m1*inv_m1m2*gmag); }
}

void InelasticCollision(Particles *particle_list, unsigned long RandomParticleType1, unsigned long RandomParticleType2, double DeltaE, bool FixedHeavyParticle)
{
	//colliding particle ititial velocities
	double v1_x = particle_list->part[RandomParticleType1].v[0];
	double v1_y = particle_list->part[RandomParticleType1].v[1];
	double v1_z = particle_list->part[RandomParticleType1].v[2];

	double v2_x = particle_list->part[RandomParticleType2].v[0];
	double v2_y = particle_list->part[RandomParticleType2].v[1];
	double v2_z = particle_list->part[RandomParticleType2].v[2];

	//particle masses
	double m1 = particle_list->part[RandomParticleType1].mass;
	double m2 = particle_list->part[RandomParticleType2].mass;
	double mr = (m1*m2)/(m1+m2);
	//center of mass velocities
	//double vcm_x = v1_x*m1+v2_x*m2;
	//double vcm_y = v1_y*m1+v2_y*m2;
	//double vcm_z = v1_z*m1+v2_z*m2;

	double du=v1_x-v2_x;
	double dv=v1_y-v2_y;
	double dw=v1_z-v2_z;

	//Cr in Bird
	double vr2 = du*du+dv*dv+dw*dw;
	double etrans = 0.5*vr2*mr;

	double divisor = 1.0 / (m1+m2);
  	double ucmf = ((m1*v1_x)+(m2*v2_x)) * divisor;
  	double vcmf = ((m1*v1_y)+(m2*v2_y)) * divisor;
  	double wcmf = ((m1*v1_z)+(m2*v2_z)) * divisor;

	double vr = sqrt(2.0 * (etrans -DeltaE)/ mr);
	if(isnan(vr))
	{ printf("# warning: negative energy encountered, etrans = %e (eV), DeltaE = %e (eV), particle type %lu with particle type %lu\n",etrans/QE,DeltaE/QE, RandomParticleType1,RandomParticleType2);
		printf("# Setting post collision v_relative to 0.0\n");
		printf("# Cross sections need to have a value of zero at and less than the threshold energy. It is recommended that you correct the cross section file.\n");
		vr=0.0;
	}

	double eps = rnd()* 2*PI;
    double cosX = 2.0*rnd() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    double ua = vr*cosX;
    double vb = vr*sinX*cos(eps);
    double wc = vr*sinX*sin(eps);

	particle_list->part[RandomParticleType1].v[0] =  ucmf + (m2*divisor)*ua;
	particle_list->part[RandomParticleType1].v[1] =  vcmf + (m2*divisor)*vb;
	particle_list->part[RandomParticleType1].v[2] =  wcmf + (m2*divisor)*wc;

	if(FixedHeavyParticle) {}
	else
	{
		particle_list->part[RandomParticleType2].v[0] = ucmf - (m1*divisor)*ua;
		particle_list->part[RandomParticleType2].v[1] = vcmf - (m1*divisor)*vb;
		particle_list->part[RandomParticleType2].v[2] = wcmf - (m1*divisor)*wc;
	}
}


void CXCollision(Particles *particle_list, unsigned long RandomParticleType1, unsigned long RandomParticleType2)
{
	//colliding particle ititial velocities
	double v1_x = particle_list->part[RandomParticleType1].v[0];
	double v1_y = particle_list->part[RandomParticleType1].v[1];
	double v1_z = particle_list->part[RandomParticleType1].v[2];

	double v2_x = particle_list->part[RandomParticleType2].v[0];
	double v2_y = particle_list->part[RandomParticleType2].v[1];
	double v2_z = particle_list->part[RandomParticleType2].v[2];

	//swap particle velocities between particles, assume that ion is the fast particle
	particle_list->part[RandomParticleType1].v[0] = v2_x;
	particle_list->part[RandomParticleType1].v[1] = v2_y;
	particle_list->part[RandomParticleType1].v[2] = v2_z;

	particle_list->part[RandomParticleType2].v[0] = v1_x;
	particle_list->part[RandomParticleType2].v[1] = v1_y;
	particle_list->part[RandomParticleType2].v[2] = v1_z;
}


double ParkElasticCollisionCosX(double eps, std::vector<double> a)
{
    if (!(a.size() > 0)) {
        printf("differential model error");
    }
    // Return cosX
    double C = a[0] + a[1]*tanh(a[2]*(eps - a[3]));
    double ef = a[4]/(eps + a[5]) + a[6];
    double eb = a[7]/(eps + a[8]) + a[9];

    double R = rnd();

    double c1 = R*ef + eb*(ef + 1 - R);
    double c2 = ef*eb + C*(ef + eb + 1) + R*(ef - eb - 1);
    double c3 = eb + R - C*(ef + eb + 1);
    double c4 = 4*R*ef*(eb + 1);

    double cosX = (c1 - sqrt(c2*c2 + c3*c4))/c3;
    return cosX;
}

double MurphyElasticCollisionCosX(double eps, std::vector<double> a)
{
    // Return cosX
    double ef = a[0]*exp(a[1]);
    double eb = a[2]*exp(a[3]);
    double N = a[4] + a[5]*tanh(a[6]*(eps-a[7]));
    double R = rnd();

    double c1 = R*ef + eb*(ef + 1 - R);
    double c2 = ef*eb + N*(ef + eb + 1) + R*(ef - eb - 1);
    double c3 = eb + R - N*(ef + eb + 1);
    double c4 = 4*R*ef*(eb + 1);

    double cosX = (c1 - sqrt(c2*c2 + c3*c4))/c3;
    return cosX;
}

// Separate the rotation and projection matrices.
void Project(double *v, double cos_chi, double sin_chi, double cos_eta, double sin_eta) {
    // Project vector into cartesian voordinates from spherical coordinates.
    v[0] = cos_chi;
    v[1] = sin_chi*cos_eta;
    v[2] = sin_chi*sin_eta;
}

void RotateFrame(double *v, double ct, double st, double cp, double sp) {
    // Compose vector *v with double-angle Euler rotation matrix
    // Copy over initial values
    double v0 = v[0];
    double v1 = v[1];
    double v2 = v[2];
    v[0] = v0*ct - v1*st;
    v[1] = v0*st*cp + v1*ct*cp - v2*sp;
    v[2] = v0*st*sp + v1*ct*sp + v2*cp;
}

void vadd_scale(double *v, double *v1, double *v2, double scale) {
    v[0] = v1[0] + scale*v2[0];
    v[1] = v1[1] + scale*v2[1];
    v[2] = v1[2] + scale*v2[2];
}
