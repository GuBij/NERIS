 #include <cmath>
 #include "particle.h"

 particle::particle (int procid, unsigned long id, double* dictParams, meteorology* wind_, meteorology* eddyDiff_, int& warmUpPeriod_, bool Langevin)
 :
  wind(wind_),
  eddyDiff(eddyDiff_),
  windParamFile("output/PARAMS_" + wind->name()+ "_" + wind->meteofile()),
  eddyDiffParamFile("output/PARAMS_" + eddyDiff->name()+ "_" + eddyDiff->meteofile()),
  windparams_ptr(&windParamFile),
  eddydiffparams_ptr(&eddyDiffParamFile),
  xMax(dictParams[4]),
  yMax(dictParams[5]),
  zMax(dictParams[6]),
  mfreq(dictParams[2]),
  dt(0),
  ppr(int(dictParams[3]/2) + (int(dictParams[3]) % 2)), // number of particles per release divided by two and rounded up
  timeIter(unsigned(floor(id/ppr))),
  releaseTime(0),
  simTime(0),
  measurementIndex(0),
  nom(wind->noRecords_()-1),
  dead(false),
  warmUpPeriod(&warmUpPeriod_),
  warmingUp(false),
  eos_(false),
  location{},
  corr{},
  up{},
  seed(1+(wind->noRecords_()+2)*procid*mfreq*ppr+id),
  generator(seed),
  useLangevinEqn(Langevin)
 {

   if (timeIter <= *warmUpPeriod)
   {
	warmingUp=true;
	eddyDiff->updateParams(eddydiffparams_ptr); wind->updateParams(windparams_ptr);
	dt = eddyDiff->horTimeStep();
	releaseTime = timeIter*dt;
   }
   else
   {
	timeIter -= *warmUpPeriod;
	int noSteps = 0;
	while ( noSteps < timeIter && measurementIndex < nom )
	{
	   eddyDiff->updateParams(eddydiffparams_ptr); wind->updateParams(windparams_ptr);
	   dt = eddyDiff->horTimeStep();
	   noSteps += mfreq/dt;
	   releaseTime += mfreq;
	   measurementIndex++;
	}
	if ( noSteps > timeIter )
	{
	   releaseTime -= (noSteps-timeIter)*dt;
	   measurementIndex--;
	}

	if (measurementIndex == nom)
	{
	  eos_ = true;
	  measurementIndex--;
	}
   }

   simTime = releaseTime;
   phi=wind->wdir(); 

   if (id % ppr == 0 && !warmingUp)
	cout << "\ntime: " << simTime << " s" << endl;

   location[2]=wind->heff();
   location[5]=location[2];

   dt_z = eddyDiff->vertTimeStep();
   nz = dt/dt_z;
   eddyDiff->corrCoeff(corr);

   Qp = (dictParams[7]*dt)/(2*ppr*dictParams[8]);
   scaleDose = dt/(2*ppr*dictParams[8]);
 }

 bool particle::inside()
 {
   if ((abs(location[0]) < xMax || abs(location[3]) < xMax) && (abs(location[1]) < yMax || abs(location[4]) < yMax) && (abs(location[2]) < zMax || abs(location[5]) < zMax))
	return true;

   if (warmingUp && (*warmUpPeriod < timeIter))
   {
	*warmUpPeriod = timeIter;
   }

   return false;
 }

 bool particle::alive()
 {
   if (dead || eos_)
	return false;

  return true;
 }

 bool particle::measurement()
 {
   if ((simTime % mfreq == 0) && (simTime > releaseTime) && !warmingUp)
     return true;

  return false; 
 }

 double particle::reflection(double zcoord, double dz, double xi, double timeStep, double elapsTime)
 {
   if ( zcoord+dz < 0 )
   {
	wind->U(Up,zcoord); eddyDiff->K(Kp,zcoord); eddyDiff->GradK(GradKp,zcoord);
	dz = (Up[2] + GradKp[2])*0.5*timeStep + sqrt(Kp[2]*timeStep)*xi + 0.25*GradKp[2]*timeStep*(xi*xi-1);
	zcoord = reflection(zcoord,dz,xi,0.5*timeStep,elapsTime);

	return zcoord;
   }
   else if ( timeStep < dt && elapsTime < dt )
   {
	zcoord += dz; elapsTime += timeStep; timeStep = dt-elapsTime;
	xi = distribution(generator);
	wind->U(Up,zcoord); eddyDiff->K(Kp,zcoord); eddyDiff->GradK(GradKp,zcoord);
        dz = (Up[2] + GradKp[2])*timeStep + sqrt(2*Kp[2]*timeStep)*xi + 0.5*GradKp[2]*timeStep*(xi*xi-1);;
	zcoord = reflection(zcoord,dz,xi,timeStep,elapsTime);

	return zcoord;
   }

   return zcoord+dz; 
 }

 void particle::eqn() 
 {
   double dW[3], dz; 
   dW[0]=distribution(generator); 
   dW[1]=distribution(generator); 
   dW[2]=distribution(generator); 

   wind->U(Up,location[2]); eddyDiff->K(Kp,location[2]); eddyDiff->GradK(GradKp,location[2]); 

   location[0] += (Up[0] + GradKp[0])*dt + sqrt(2*Kp[0]*dt)*dW[0];
   location[1] += (Up[1] + GradKp[1])*dt + sqrt(2*Kp[1]*dt)*dW[1];
   location[2] = abs(location[2] + (Up[2] + GradKp[2])*dt_z + sqrt(2*Kp[2]*dt_z)*dW[2]); 

   wind->U(Up,location[5]); eddyDiff->K(Kp,location[5]); eddyDiff->GradK(GradKp,location[5]); 

   location[3] += (Up[0] + GradKp[0])*dt - sqrt(2*Kp[0]*dt)*dW[0];
   location[4] += (Up[1] + GradKp[1])*dt - sqrt(2*Kp[1]*dt)*dW[1];
   location[5] = abs(location[5] + (Up[2] + GradKp[2])*dt_z - sqrt(2*Kp[2]*dt_z)*dW[2]); 

   unsigned i = 1;
   while (i < nz)
   {
        dW[2]=distribution(generator);
        wind->U(Up,location[2]); eddyDiff->K(Kp,location[2]); eddyDiff->GradK(GradKp,location[2]);
        location[2] = abs(location[2] + (Up[2] + GradKp[2])*dt_z + sqrt(2*Kp[2]*dt_z)*dW[2]);

        wind->U(Up,location[5]); eddyDiff->K(Kp,location[5]); eddyDiff->GradK(GradKp,location[5]);
        location[5] = abs(location[5] + (Up[2] + GradKp[2])*dt_z - sqrt(2*Kp[2]*dt_z)*dW[2]);

        i++;
   }
 }

 void particle::LangevinEqn()
 {
   double dW[3];                //brownian increments
   dW[0]=distribution(generator);
   dW[1]=distribution(generator);
   dW[2]=distribution(generator);

   wind->U(Up,location[2]); eddyDiff->K(Kp,location[2]); GradSig2_ = eddyDiff->GradSig2(location[2]);

   up[0] = corr[0]*up[0] + Kp[0]*dW[0];
   up[1] = corr[1]*up[1] + Kp[1]*dW[1];
   up[2] = corr[2]*up[2] + Kp[2]*dW[2] + dt_z*GradSig2_;

   location[0] += (Up[0] + up[0]*cos(phi)-up[1]*sin(phi))*dt;
   location[1] += (Up[1] + up[0]*sin(phi)+up[1]*cos(phi))*dt;
   location[2] = location[2] + (Up[2] + up[2])*dt_z;

   if (location[2] < 0)
   {
	location[2] *= -1; up[2] *= -1;
   }

   wind->U(Up,location[5]); eddyDiff->K(Kp,location[5]); GradSig2_ = eddyDiff->GradSig2(location[5]);

   up[3] = corr[0]*up[3] - Kp[0]*dW[0];
   up[4] = corr[1]*up[4] - Kp[1]*dW[1];
   up[5] = corr[2]*up[5] - Kp[2]*dW[2] + dt_z*GradSig2_;

   location[3] += (Up[0] + up[3]*cos(phi)-up[4]*sin(phi))*dt;
   location[4] += (Up[1] + up[3]*sin(phi)+up[4]*cos(phi))*dt;
   location[5] = location[5] + (Up[2] + up[5])*dt_z;

   if (location[5] < 0)
   {
        location[5] *= -1; up[5] *= -1;
   }

   unsigned i = 1;
   while (i < nz)
   {
	dW[2]=distribution(generator);
	wind->U(Up,location[2]); eddyDiff->K(Kp,location[2]); GradSig2_ = eddyDiff->GradSig2(location[2]);
	up[2] = corr[2]*up[2] + Kp[2]*dW[2] + dt_z*GradSig2_;
	location[2] = location[2] + (Up[2] + up[2])*dt_z;

 	if (location[2] < 0)
	{
	  location[2] *= -1; up[2] *= -1;
	}

	wind->U(Up,location[5]); eddyDiff->K(Kp,location[5]); GradSig2_ = eddyDiff->GradSig2(location[5]);
	up[5] = corr[2]*up[5] - Kp[2]*dW[2] + dt_z*GradSig2_;
	location[5] = location[5] + (Up[2] + up[5])*dt_z;

        if (location[5] < 0)
        {
          location[5] *= -1; up[5] *= -1;
        }

	i++;
   }
 }

 void particle::advance () 
 {
   if (measurement() && measurementIndex < nom-1)
   {
	measurementIndex++; 
	wind->updateParams(windparams_ptr); eddyDiff->updateParams(eddydiffparams_ptr); eddyDiff->corrCoeff(corr); phi = wind->wdir();
	dt = eddyDiff->horTimeStep(); dt_z = eddyDiff->vertTimeStep(); nz = dt/dt_z;
   }
   else if (measurement())
  	 dead=true; 

   if (useLangevinEqn)
	LangevinEqn();
   else
	eqn(); 
 
   timeIter++; simTime += dt;

   if ((releaseTime > 0) && warmingUp && (timeIter > *warmUpPeriod))
   {
	warmingUp = false; releaseTime = timeIter*dt - releaseTime; timeIter -= *warmUpPeriod; simTime = timeIter*dt; releaseTime = timeIter*dt - releaseTime; 
   }
 }

