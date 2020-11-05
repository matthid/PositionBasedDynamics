#include "ExtendedPositionBasedDynamics.h"
#include "MathFunctions.h"
#include <cfloat>

using namespace PBD;

const Real PBD::ExtendedPositionBasedDynamics::eps = static_cast<Real>(1e-6);
const Real PBD::ExtendedPositionBasedDynamics::maxRotationPerSubstep = 0.5;

//////////////////////////////////////////////////////////////////////////
// ExtendedPositionBasedDynamics
//////////////////////////////////////////////////////////////////////////
const Vector3r ExtendedPositionBasedDynamics::SOFTBODY_INVERSE_INERTIA = Vector3r(0,0,0);


bool PBD::ExtendedPositionBasedDynamics::calc_lagrange_PositionalConstraint(
	const Vector3r& p_corr, 
	const Vector3r& p0, const Real invMass0, const Vector3r& inverseInertia0, const Quaternionr& inverseRotation0,
	const Vector3r& p1, const Real invMass1, const Vector3r& inverseInertia1, const Quaternionr& inverseRotation1,
	const Real compliance, const Real dt,
	Real& lagrangeMultiplier,
	Vector3r& positionalImpulse,
	Vector3r& positionalImpulse0,
	Vector3r& positionalImpulse1)
{
	Real magnitude = p_corr.norm();
	if (abs(magnitude) < eps) {
		return false;
	}
	Vector3r direction = p_corr.normalized();
	Vector3r direction0 = inverseRotation0 * direction;
	Vector3r direction1 = inverseRotation1 * direction;

	Vector3r p0_cross = p0.cross(direction0);
	Real generalizedMass0 = invMass0 + (p0_cross.transpose() * inverseInertia0.asDiagonal()).dot(p0_cross);

	Vector3r p1_cross = p1.cross(direction1);
	Real generalizedMass1 = invMass1 + (p1_cross.transpose() * inverseInertia1.asDiagonal()).dot(p1_cross);

	Real w = generalizedMass0 + generalizedMass1;
	if (abs(w) < eps) {
		return false;
	}

	Real alpha = compliance / (dt * dt);
	Real lagrangeMultiplierCorr = (-magnitude - alpha * lagrangeMultiplier) / (w + alpha);
	
	lagrangeMultiplier += lagrangeMultiplierCorr;
	positionalImpulse = lagrangeMultiplierCorr * direction;
	positionalImpulse0 = lagrangeMultiplierCorr * direction0;
	positionalImpulse1 = lagrangeMultiplierCorr * direction1;
	return true;
}

bool PBD::ExtendedPositionBasedDynamics::update_positions(
	Vector3r& p0, const Real invMass0, 
	Vector3r& p1, const Real invMass1,
	const Vector3r positionalImpuse)
{
	p0 += invMass0 * positionalImpuse;
	p1 -= invMass1 * positionalImpuse;
	return false;
}

bool PBD::ExtendedPositionBasedDynamics::scale_rotation_update(
	const Vector3r inverseInertia0,
	const Vector3r inverseInertia1,
	const Vector3r rotationUpdate0,
	const Vector3r rotationUpdate1,
	Vector3r& scaledRotationUpdate0,
	Vector3r& scaledRotationUpdate1)
{
	Vector3r rot_update0 = inverseInertia0.asDiagonal() * rotationUpdate0;
	scaledRotationUpdate0 = 0.5 * rot_update0;

	Vector3r rot_update1 = inverseInertia1.asDiagonal() * rotationUpdate1;
	scaledRotationUpdate1 = 0.5 * rot_update1;
	return true;
}

bool PBD::ExtendedPositionBasedDynamics::update_velocities(
	Quaternionr& orientiation0, 
	Quaternionr& orientiation1,
	const Vector3r scaledRotationUpdate0,
	const Vector3r scaledRotationUpdate1)
{
	Quaternionr angVel0(0.0, scaledRotationUpdate0[0], scaledRotationUpdate0[1], scaledRotationUpdate0[2]);
	orientiation0.coeffs() += (angVel0 * orientiation0).coeffs();
	orientiation0.normalize();

	Quaternionr angVel1(0.0, scaledRotationUpdate1[0], scaledRotationUpdate1[1], scaledRotationUpdate1[2]);
	orientiation1.coeffs() -= (angVel1 * orientiation1).coeffs();
	orientiation1.normalize();
	return true;
}

bool PBD::ExtendedPositionBasedDynamics::calc_lagrange_AngularConstraint(
	const Vector3r& q_corr,
	const Vector3r& inverseInertia0, const Quaternionr& inverseRotation0,
	const Vector3r& inverseInertia1, const Quaternionr& inverseRotation1,
	const Real compliance, const Real dt, 
	Real& lagrangeMultiplier,
	Vector3r& positionalImpulse0,
	Vector3r& positionalImpulse1)
{
	Real magnitude = q_corr.norm(); // rotation angle
	if (abs(magnitude) < eps) {
		return false;
	}

	Vector3r direction = q_corr.normalized();
	Vector3r direction0 = inverseRotation0 * direction;
	Vector3r direction1 = inverseRotation1 * direction;

	Real generalizedMass0 = (direction0.transpose() * inverseInertia0.asDiagonal()).dot(direction0);
	Real generalizedMass1 = (direction1.transpose() * inverseInertia1.asDiagonal()).dot(direction1);

	Real w = generalizedMass0 + generalizedMass1;
	if (abs(w) < eps) {
		return false;
	}
	Real alpha = compliance / (dt * dt);
	Real lagrangeMultiplierCorr = (-magnitude - alpha * lagrangeMultiplier) / (w + alpha);

	lagrangeMultiplier += lagrangeMultiplierCorr;
	positionalImpulse0 = lagrangeMultiplierCorr * direction0;
	positionalImpulse1 = lagrangeMultiplierCorr * direction1;
	return true;
}

bool PBD::ExtendedPositionBasedDynamics::limitAngle(
	//int body0, int body1, int n, int a, int b, int minAngle, int maxAngle, int compliance, Real dt, Real maxCorr)
	const Vector3r& common_rot_axis, const Vector3r& axis0, const Vector3r& axis1, const Real minAngle, const Real maxAngle, Vector3r& apply_rot)
{

	Real phi = asin(axis0.cross(axis1).dot(common_rot_axis));
	if (axis0.dot(axis1) < 0) {
		phi = 2 * M_PI - phi;
	}

	if (phi > M_PI) {
		phi -= 2 * M_PI;
	}

	if (phi < -M_PI) {
		phi += 2 * M_PI;
	}

	if (phi < minAngle || phi > maxAngle) {
		//phi = clamp(phi, minAngle, maxAngle);
		phi = std::max(minAngle, std::min(phi, maxAngle));
		Quaternionr rot(AngleAxisr(phi, common_rot_axis));
		Vector3r rot_axis0 = rot * axis0;
		apply_rot = rot_axis0.cross(axis1);
		return true;
	}

	return false;
}

