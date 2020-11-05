#ifndef EXTENDED_POSITION_BASED_DYNAMICS_H
#define EXTENDED_POSITION_BASED_DYNAMICS_H

#include "Common/Common.h"
#define _USE_MATH_DEFINES
#include <math.h>

// ------------------------------------------------------------------------------------
namespace PBD
{
	class ExtendedPositionBasedDynamics
	{
	public:
		static const Real eps;

		// -------------- XPBD -----------------------------------------------------
		static const Vector3r SOFTBODY_INVERSE_INERTIA;

		template<class TBody>
		static Vector3r& getBodyVel(TBody& obj)
		{
			return obj.getVelocity();
		}

		template<class TBody>
		static Vector3r& getBodyPos(TBody& obj)
		{
			return obj.getPosition();
		}

		template<class TBody>
		static Vector3r& getBodyAccel(TBody& obj)
		{
			return obj.getAcceleration();
		}

		template<class TBody>
		static Quaternionr* getBodyRot(TBody& obj)
		{
			return &obj.getRotation();
		}

		//template<class TBody>
		//static Quaternionr* getBodyRot(...)
		//{
		//	return nullptr;
		//}

		template<class TBody>
		static const Vector3r* getBodyInertiaInv(TBody& obj)
		{
			return &obj.getInertiaTensorInverse();
		}

		//template<class TBody>
		//static Vector3r* getBodyInertiaInv(...)
		//{
		//	return nullptr;
		//}

		template<class TBody>
		static Real getBodyMassInv(TBody& obj)
		{
			return obj.getInvMass();
		}

		//template<class TBody>
		//static Real getBodyMassInv(...)
		//{
		//	return 0;
		//}

		template <class TBody>
		static Real getInverseMass(TBody& body, const Vector3r& normal, const Vector3r* pos = nullptr) {
			Vector3r n;
			if (pos == nullptr) {
				n = normal;
			}
			else {
				Vector3r diff = *pos - getBodyPos(body);
				n = diff.cross(normal);
			}

			Quaternionr* rot = getBodyRot(body);
			n = rot->conjugate() * n;

			Real w = 0;
			const Vector3r* invInertia = getBodyInertiaInv(body);
			if (invInertia != nullptr) {
				w = (normal.transpose() * invInertia->asDiagonal()).dot(normal);
			}
			if (pos != nullptr) {
				w += getBodyMassInv(body);
			}

			return w;
		}

		static const Real maxRotationPerSubstep;

		template <class TBody>
		static void applyRotation(TBody& body, const Vector3r& rot, Real scale = 1) {
			Real maxPhi = 0.5;
			Real phi = rot.norm();
			if (phi * scale > maxRotationPerSubstep) {
				scale = maxRotationPerSubstep / phi;
			}

			Quaternionr* bodyRot = getBodyRot(body);
			if (bodyRot != nullptr) {
				Quaternionr dq = Quaternionr(0, rot.x() * scale, rot.y() * scale, rot.z() * scale) * *bodyRot;
				bodyRot->x() += 0.5 * dq.x();
				bodyRot->y() += 0.5 * dq.y();
				bodyRot->z() += 0.5 * dq.z();
				bodyRot->w() += 0.5 * dq.w();
				bodyRot->normalize();
			}
		}

		template <class TBody>
		static void applyCorrection(TBody& body, const Vector3r& corr, const Vector3r* pos = nullptr, const bool velocityLevel = false) {
			Vector3r dq;
			if (pos == nullptr) {
				dq = corr;
			}
			else {
				Vector3r scaled = getBodyMassInv(body) * corr;
				if (velocityLevel) {
					getBodyVel(body) += scaled;
				}
				else {
					getBodyPos(body) += scaled;
				}
				Vector3r diff = *pos - getBodyPos(body);
				dq = diff.cross(corr);
			}

			Quaternionr* rot = getBodyRot(body);
			dq = rot->conjugate() * dq;
			Vector3r invInertia = *getBodyInertiaInv(body);
			dq = *rot * Vector3r(invInertia.x() * dq.x(), invInertia.y() * dq.y(), invInertia.z() * dq.z());
			if (velocityLevel) {
				getBodyAccel(body) += dq;
			}
			else {
				applyRotation(body, dq);
			}
		}


		template <class TBody>
		static void applyBodyPairCorrection(
				TBody& body0, TBody& body1, const Vector3r& corr, const Real compliance, const Real dt, const Vector3r* pos0 = nullptr, const Vector3r* pos1 = nullptr, const bool velocityLevel = false)
		{
			Real length = corr.norm();
			if (abs(length) < eps) {
				return;
			}

			Vector3r normal = corr;
			normal.normalize();
			Real w0 = getInverseMass(body0, normal, pos0);
			Real w1 = getInverseMass(body1, normal, pos1);
			Real w = w0 + w1;
			if (abs(w) < eps) {
				return;
			}

			Real lambda = -length / (w + compliance / dt / dt);
			normal *= -lambda;
			applyCorrection(body0, normal, pos0, velocityLevel);
			normal *= -1;
			applyCorrection(body1, normal, pos1, velocityLevel);
		}


		/** See 3.3.1 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Formula (2), (3), (4), (5)
		*/
		static bool calc_lagrange_PositionalConstraint(
			const Vector3r& p_corr,
			const Vector3r& p0, const Real invMass0, const Vector3r& inverseInertia0, const Quaternionr& inverseRotation0,
			const Vector3r& p1, const Real invMass1, const Vector3r& inverseInertia1, const Quaternionr& inverseRotation1,
			const Real compliance,  // inverse stiffness
			const Real dt,
			Real& lagrangeMultiplier,
			Vector3r& positionalImpulse,
			Vector3r& positionalImpulse0,
			Vector3r& positionalImpulse1);

		/** See 3.3.1 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Formula (6) & (7)
		*/
		static bool update_positions(
			Vector3r& p0, const Real invMass0,
			Vector3r& p1, const Real invMass1,
			const Vector3r positionalImpuse);

		/** See 3.3.1 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Formula (8) & (9) -> Use ( pos_n x positionalImpuse_n ) as rotationUpdate for positional updates
		* Formula (15) & (16) -> Use positionalUpdate as rotiationUpdate for angular velocity updates
		*/
		static bool scale_rotation_update(
			const Vector3r inverseInertia0,
			const Vector3r inverseInertia1,
			const Vector3r rotationUpdate0,
			const Vector3r rotationUpdate1,
			Vector3r& scaledRotationUpdate0,
			Vector3r& scaledRotationUpdate1);

		/** See 3.3.1 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Use output of scale_rotation_update
		*/
		static bool update_velocities(
			Quaternionr& orientiation0,
			Quaternionr& orientiation1,
			const Vector3r scaledRotationUpdate0,
			const Vector3r scaledRotationUpdate1);

		/** See 3.3.2 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Formula (11), (12), (13), (14)
		 */
		static bool calc_lagrange_AngularConstraint(
			const Vector3r& q_corr,
			const Vector3r& inverseInertia0, // rest space/diagonal
			const Quaternionr& inverseRotation0,
			const Vector3r& inverseInertia1, // rest space/diagonal
			const Quaternionr& inverseRotation1,
			const Real compliance, // inverse stiffness
			const Real dt,
			Real& lagrangeMultiplier,
			Vector3r& positionalImpulse0,
			Vector3r& positionalImpulse1);

	private:
		// https://github.com/matthias-research/pages/blob/fb41d362b229189a55143e2968e844a7cf5464c6/challenges/PBD.js#L219-L250
		/** See 3.4 in Detailed Rigid Body Simulation with Extended Position Based Dynamics (2020)
		* Algorithm 3
		 */
		static bool limitAngle(
			const Vector3r& common_rot_axis, const Vector3r& axis0, const Vector3r& axis1, const Real minAngle, const Real maxAngle, Vector3r& apply_rot
		);

	};
}

#endif
