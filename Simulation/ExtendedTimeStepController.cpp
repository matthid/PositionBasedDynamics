#include "ExtendedTimeStepController.h"
#include "Simulation/TimeManager.h"
#include "PositionBasedDynamics/PositionBasedRigidBodyDynamics.h"
#include "PositionBasedDynamics/TimeIntegration.h"
#include <iostream>
#include "PositionBasedDynamics/PositionBasedDynamics.h"
#include "Utils/Timing.h"

using namespace PBD;
using namespace std;
using namespace GenParam;

// int ExtendedTimeStepController::SOLVER_ITERATIONS = -1;
// int ExtendedTimeStepController::SOLVER_ITERATIONS_V = -1;
int ExtendedTimeStepController::NUM_SUBSTEPS = -1;
//int ExtendedTimeStepController::MAX_ITERATIONS = -1;
//int ExtendedTimeStepController::MAX_ITERATIONS_V = -1;
int ExtendedTimeStepController::VELOCITY_UPDATE_METHOD = -1;
int ExtendedTimeStepController::ENUM_VUPDATE_FIRST_ORDER = -1;
int ExtendedTimeStepController::ENUM_VUPDATE_SECOND_ORDER = -1;


ExtendedTimeStepController::ExtendedTimeStepController()
{
	m_velocityUpdateMethod = 0;
	m_iterations = 0;
	m_iterationsV = 0;
	m_numSubsteps = 1;
	m_maxIterations = 1;
	m_maxIterationsV = 1;
	m_collisionDetection = NULL;
}

ExtendedTimeStepController::~ExtendedTimeStepController(void)
{
}

void ExtendedTimeStepController::initParameters()
{
	TimeStep::initParameters();

	// 	SOLVER_ITERATIONS = createNumericParameter("iterations", "Iterations", &m_iterations);
	// 	setGroup(SOLVER_ITERATIONS, "PBD");
	// 	setDescription(SOLVER_ITERATIONS, "Iterations required by the solver.");
	// 	getParameter(SOLVER_ITERATIONS)->setReadOnly(true);

	NUM_SUBSTEPS = createNumericParameter("numSubsteps", "Number of substeps", &m_numSubsteps);
	setGroup(NUM_SUBSTEPS, "XPBD");
	setDescription(NUM_SUBSTEPS, "Number of substeps to use.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(NUM_SUBSTEPS))->setMinValue(1);

	//MAX_ITERATIONS = createNumericParameter("maxIterations", "Max. iterations", &m_maxIterations);
	//setGroup(MAX_ITERATIONS, "XPBD");
	//setDescription(MAX_ITERATIONS, "Maximal number of iterations of the solver.");
	//static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	// 	SOLVER_ITERATIONS_V = createNumericParameter("iterationsV", "Velocity iterations", &m_iterationsV);
	// 	setGroup(SOLVER_ITERATIONS_V, "PBD");
	// 	setDescription(SOLVER_ITERATIONS_V, "Iterations required by the velocity solver.");
	// 	getParameter(SOLVER_ITERATIONS_V)->setReadOnly(true);

	//MAX_ITERATIONS_V = createNumericParameter("maxIterationsV", "Max. velocity iterations", &m_maxIterationsV);
	//setGroup(MAX_ITERATIONS_V, "XPBD");
	//setDescription(MAX_ITERATIONS_V, "Maximal number of iterations of the velocity solver.");
	//static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(0);

	VELOCITY_UPDATE_METHOD = createEnumParameter("velocityUpdateMethod", "Velocity update method", &m_velocityUpdateMethod);
	setGroup(VELOCITY_UPDATE_METHOD, "XPBD");
	setDescription(VELOCITY_UPDATE_METHOD, "Velocity method.");
	EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(VELOCITY_UPDATE_METHOD));
	enumParam->addEnumValue("First Order Update", ENUM_VUPDATE_FIRST_ORDER);
	enumParam->addEnumValue("Second Order Update", ENUM_VUPDATE_SECOND_ORDER);
}

void ExtendedTimeStepController::step(SimulationModel& model)
{
	START_TIMING("simulation step");
	TimeManager* tm = TimeManager::getCurrent();
	const Real timestep = tm->getTimeStepSize();

	// See Detailed Rigid Body Simulation with Extended Position Based Dynamics, Algorithm 2
	if (m_collisionDetection)
	{
		START_TIMING("collision detection");
		m_collisionDetection->collisionDetection(model);
		STOP_TIMING_AVG;
	}

	const Real h = timestep / m_numSubsteps;

	// I assume this is for changes in gravity?
	clearAccelerations(model);

	SimulationModel::RigidBodyVector& rb = model.getRigidBodies();
	ParticleData& pd = model.getParticles();
	OrientationData& od = model.getOrientations();

	const int numBodies = (int)rb.size();

	for (int s = 0; s < m_numSubsteps; s++)
	{
		LOG_CALC(LOG_INFO << "Simuation Step (Substep: " << s << "): " << timestep << ", devided in " << h);
#pragma omp parallel if(numBodies > MIN_PARALLEL_SIZE) default(shared)
		{
			//////////////////////////////////////////////////////////////////////////
			// rigid body model
			//////////////////////////////////////////////////////////////////////////
#pragma omp for schedule(static) nowait
			for (int i = 0; i < numBodies; i++)
			{
				// x_prev = x
				rb[i]->getLastPosition() = rb[i]->getOldPosition();
				rb[i]->getOldPosition() = rb[i]->getPosition();
				// v <- v + h*f_ext/m
				// x <- x + hv
				TimeIntegration::semiImplicitEuler(h, rb[i]->getMass(), rb[i]->getPosition(), rb[i]->getVelocity(), rb[i]->getAcceleration());
				
				// q_prev = q
				rb[i]->getLastRotation() = rb[i]->getOldRotation();
				rb[i]->getOldRotation() = rb[i]->getRotation();
				// w <- w + h * I^-1 * 
				// q <- q + h * 0.5 * [w_x, w_y, w_z, 0] * q
				// normalize: q <- q / |q|
				TimeIntegration::semiImplicitEulerRotation(h, rb[i]->getMass(), rb[i]->getInertiaTensorInverseW(), rb[i]->getRotation(), rb[i]->getAngularVelocity(), rb[i]->getTorque());
				//rb[i]->rotationUpdated();

				LOG_CALC(LOG_INFO << "Updated body: " << rb[i]->getName()
					<< "\n\tx_prev: " << rb[i]->getOldPosition() << ", x: " << rb[i]->getPosition() 
					<< "\n\tv: " << rb[i]->getVelocity() << ", w: " << rb[i]->getAngularVelocity()
					<< "\n\tq_prev: " << rb[i]->getOldRotation() << ", q: " << rb[i]->getRotation());
			}

			//////////////////////////////////////////////////////////////////////////
			// particle model
			// x_prev = x
			// v <- v + h*f_ext/m
			// x <- x + hv
			//////////////////////////////////////////////////////////////////////////
#pragma omp for schedule(static) 
			for (int i = 0; i < (int)pd.size(); i++)
			{
				// x_prev = x
				pd.getLastPosition(i) = pd.getOldPosition(i);
				pd.getOldPosition(i) = pd.getPosition(i);
				// v <- v + h*f_ext/m
				// x <- x + hv
				TimeIntegration::semiImplicitEuler(h, pd.getMass(i), pd.getPosition(i), pd.getVelocity(i), pd.getAcceleration(i));
			}
		}

		START_TIMING("position constraints projection");
		// SolvePositions
		positionConstraintProjection(model, h);
		STOP_TIMING_AVG;

#pragma omp parallel if(numBodies > MIN_PARALLEL_SIZE) default(shared)
		{
			// Update velocities for rigid bodies
#pragma omp for schedule(static) nowait
			for (int i = 0; i < numBodies; i++)
			{
				if (m_velocityUpdateMethod == 0)
				{
					// v <- x - x_prev / h
					TimeIntegration::velocityUpdateFirstOrder(h, rb[i]->getMass(), rb[i]->getPosition(), rb[i]->getOldPosition(), rb[i]->getVelocity());

					// dq <- q * q_prev^-1
					// w <- 2 * [d q_x, d q_y, d q_z] / h
					// w <- (d q_w >= 0) ? w : -w
					TimeIntegration::angularVelocityUpdateFirstOrder(h, rb[i]->getMass(), rb[i]->getRotation(), rb[i]->getOldRotation(), rb[i]->getAngularVelocity(), true);
				}
				else
				{
					TimeIntegration::velocityUpdateSecondOrder(h, rb[i]->getMass(), rb[i]->getPosition(), rb[i]->getOldPosition(), rb[i]->getLastPosition(), rb[i]->getVelocity());
					TimeIntegration::angularVelocityUpdateSecondOrder(h, rb[i]->getMass(), rb[i]->getRotation(), rb[i]->getOldRotation(), rb[i]->getLastRotation(), rb[i]->getAngularVelocity());
				}


				LOG_CALC(LOG_INFO << "Updated body after positional solve: " << rb[i]->getName()
					<< "\n\tv: " << rb[i]->getVelocity()
					<< "\n\tx_prev: " << rb[i]->getOldPosition() << ", x: " << rb[i]->getPosition()
					<< "\n\tw: " << rb[i]->getAngularVelocity()
					<< "\n\tq_prev: " << rb[i]->getOldRotation() << ", q: " << rb[i]->getRotation());
			}

			// Update velocities	
#pragma omp for schedule(static) 
			for (int i = 0; i < (int)pd.size(); i++)
			{
				if (m_velocityUpdateMethod == 0)
					// v <- x - x_prev / h
					TimeIntegration::velocityUpdateFirstOrder(h, pd.getMass(i), pd.getPosition(i), pd.getOldPosition(i), pd.getVelocity(i));
				else
					TimeIntegration::velocityUpdateSecondOrder(h, pd.getMass(i), pd.getPosition(i), pd.getOldPosition(i), pd.getLastPosition(i), pd.getVelocity(i));
			}

			// Update velocites of orientations
#pragma omp for schedule(static) 
			for (int i = 0; i < (int)od.size(); i++)
			{
				if (m_velocityUpdateMethod == 0)
					// dq <- q * q_prev^-1
					// w <- 2 * [d q_x, d q_y, d q_z] / h
					// w <- (d q_w >= 0) ? w : -w
					TimeIntegration::angularVelocityUpdateFirstOrder(h, od.getMass(i), od.getQuaternion(i), od.getOldQuaternion(i), od.getVelocity(i), true);
				else
					TimeIntegration::angularVelocityUpdateSecondOrder(h, od.getMass(i), od.getQuaternion(i), od.getOldQuaternion(i), od.getLastQuaternion(i), od.getVelocity(i));
			}
		}


		START_TIMING("velocity constraints projection");
		// SolveVelocities
		velocityConstraintProjection(model);
		STOP_TIMING_AVG;


		// Update Inertia Matrices for all bodies (we do this once per substep only)
		// We only need this at the beginning `semiImplicitEulerRotation` and use rest state when solving constraits
#pragma omp for schedule(static) nowait
		for (int i = 0; i < numBodies; i++)
		{
			rb[i]->rotationUpdated();
		}
	}

#pragma omp parallel if(numBodies > MIN_PARALLEL_SIZE) default(shared)
#pragma omp for schedule(static) nowait
	// Update Geometry?
	for (int i = 0; i < numBodies; i++)
	{
		// update geometry
		if (rb[i]->getMass() != 0.0)
			rb[i]->getGeometry().updateMeshTransformation(rb[i]->getPosition(), rb[i]->getRotationMatrix());
	}

	//////////////////////////////////////////////////////////////////////////
	// update motor joint targets
	//////////////////////////////////////////////////////////////////////////
	SimulationModel::ConstraintVector& constraints = model.getConstraints();
	for (unsigned int i = 0; i < constraints.size(); i++)
	{
		if ((constraints[i]->getTypeId() == TargetAngleMotorHingeJoint::TYPE_ID) ||
			(constraints[i]->getTypeId() == TargetVelocityMotorHingeJoint::TYPE_ID) ||
			(constraints[i]->getTypeId() == TargetPositionMotorSliderJoint::TYPE_ID) ||
			(constraints[i]->getTypeId() == TargetVelocityMotorSliderJoint::TYPE_ID))
		{
			MotorJoint* motor = (MotorJoint*)constraints[i];
			const std::vector<Real> sequence = motor->getTargetSequence();
			if (sequence.size() > 0)
			{
				Real time = tm->getTime();
				const Real sequenceDuration = sequence[sequence.size() - 2] - sequence[0];
				if (motor->getRepeatSequence())
				{
					while (time > sequenceDuration)
						time -= sequenceDuration;
				}
				unsigned int index = 0;
				while ((2 * index < sequence.size()) && (sequence[2 * index] <= time))
					index++;

				// linear interpolation
				Real target = 0.0;
				if (2 * index < sequence.size())
				{
					const Real alpha = (time - sequence[2 * (index - 1)]) / (sequence[2 * index] - sequence[2 * (index - 1)]);
					target = (static_cast<Real>(1.0) - alpha) * sequence[2 * index - 1] + alpha * sequence[2 * index + 1];
				}
				else
					target = sequence[sequence.size() - 1];
				motor->setTarget(target);
			}
		}
	}

	// compute new time	
	tm->setTime(tm->getTime() + timestep);
	STOP_TIMING_AVG;
}

void ExtendedTimeStepController::reset()
{
	m_iterations = 0;
	m_iterationsV = 0;
	m_numSubsteps = 5;
	m_maxIterations = 1;
	m_maxIterationsV = 1;
}

void ExtendedTimeStepController::positionConstraintProjection(SimulationModel& model, const Real dt)
{
	m_iterations = 0;

	// init constraint groups if necessary
	model.initConstraintGroups();

	SimulationModel::RigidBodyVector& rb = model.getRigidBodies();
	SimulationModel::ConstraintVector& constraints = model.getConstraints();
	SimulationModel::ConstraintGroupVector& groups = model.getConstraintGroups();
	SimulationModel::RigidBodyContactConstraintVector& contacts = model.getRigidBodyContactConstraints();
	SimulationModel::ParticleSolidContactConstraintVector& particleTetContacts = model.getParticleSolidContactConstraints();

	// init constraints for this time step if necessary
	for (auto& constraint : constraints)
	{
		constraint->initConstraintBeforeProjection(model);
	}


	//while (m_iterations < m_maxIterations)
	//{
		for (unsigned int group = 0; group < groups.size(); group++)
		{
			const int groupSize = (int)groups[group].size();
#pragma omp parallel if(groupSize > MIN_PARALLEL_SIZE) default(shared)
			{
#pragma omp for schedule(static) 
				for (int i = 0; i < groupSize; i++)
				{
					const unsigned int constraintIndex = groups[group][i];

					//constraints[constraintIndex]->updateConstraint(model);
					constraints[constraintIndex]->extendedPBDRigidBodyUpdate(model, dt);
					//constraints[constraintIndex]->solvePositionConstraint(model, m_iterations);
				}
			}
		}

		for (unsigned int i = 0; i < particleTetContacts.size(); i++)
		{
			particleTetContacts[i].extendedPBDRigidBodyUpdate(model, dt);
		}

	//	m_iterations++;
	//}
}


void ExtendedTimeStepController::velocityConstraintProjection(SimulationModel& model)
{

	m_iterationsV = 0;

	// init constraint groups if necessary
	model.initConstraintGroups();

	SimulationModel::RigidBodyVector& rb = model.getRigidBodies();
	SimulationModel::ConstraintVector& constraints = model.getConstraints();
	SimulationModel::ConstraintGroupVector& groups = model.getConstraintGroups();
	SimulationModel::RigidBodyContactConstraintVector& rigidBodyContacts = model.getRigidBodyContactConstraints();
	SimulationModel::ParticleRigidBodyContactConstraintVector& particleRigidBodyContacts = model.getParticleRigidBodyContactConstraints();
	SimulationModel::ParticleSolidContactConstraintVector& particleTetContacts = model.getParticleSolidContactConstraints();

	// TODO:
	// for each contact do
	//   (29) v <- (v_1 + omega_1 + r_1) - (v_2 + omega_2 + r_2)
	//   (29) v_n <- dot(n, v)
	//   (29) v_t <- v - n*v_n
	//   (30) dv <- -(v_t / | v_t |)*min(h * mü_d * f_n, |v_t|)
	//        (mü_d = dynamic friction coefficient, f_n = lambda_n / dt^2 (normal force)
	//   (31) dv <- (v2 - v1)*min(mü_lin * dt, 1)
	//   (32) domega <- (omega_2 - omega_1)*min(mü_ang * dt, 1)
	//   (32) domega <- (omega_2 - omega_1)*min(mü_ang * dt, 1)
	//   (33)


	for (unsigned int group = 0; group < groups.size(); group++)
	{
		const int groupSize = (int)groups[group].size();
#pragma omp parallel if(groupSize > MIN_PARALLEL_SIZE) default(shared)
		{
#pragma omp for schedule(static) 
			for (int i = 0; i < groupSize; i++)
			{
				const unsigned int constraintIndex = groups[group][i];
				//constraints[constraintIndex]->updateConstraint(model);
			}
		}
	}

	//while (m_iterationsV < m_maxIterationsV)
	//{
		for (unsigned int group = 0; group < groups.size(); group++)
		{
			const int groupSize = (int)groups[group].size();
#pragma omp parallel if(groupSize > MIN_PARALLEL_SIZE) default(shared)
			{
#pragma omp for schedule(static) 
				for (int i = 0; i < groupSize; i++)
				{
					const unsigned int constraintIndex = groups[group][i];
					//constraints[constraintIndex]->solveVelocityConstraint(model, m_iterationsV);
				}
			}
		}

		// solve contacts
		for (unsigned int i = 0; i < rigidBodyContacts.size(); i++)
		{
		//	rigidBodyContacts[i].solveVelocityConstraint(model, m_iterationsV);
		}
		for (unsigned int i = 0; i < particleRigidBodyContacts.size(); i++)
		{
		//	particleRigidBodyContacts[i].solveVelocityConstraint(model, m_iterationsV);
		}
		for (unsigned int i = 0; i < particleTetContacts.size(); i++)
		{
		//	particleTetContacts[i].solveVelocityConstraint(model, m_iterationsV);
		}
		m_iterationsV++;
	//}
}


