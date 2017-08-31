#ifndef __PARTICLEDATA_H__
#define __PARTICLEDATA_H__

#include <vector>
#include "Common/Common.h"


namespace PBD
{
	/** This class encapsulates the state of all vertices.
	* All parameters are stored in individual arrays.
	*/
	class VertexData
	{
	private:
		std::vector<Vector3r> m_x;

	public:
		FORCE_INLINE VertexData(void) :
			m_x()
		{
		}

		FORCE_INLINE ~VertexData(void)
		{
			m_x.clear();
		}

		FORCE_INLINE void addVertex(const Vector3r &vertex)
		{
			m_x.push_back(vertex);
		}

		FORCE_INLINE Vector3r &getPosition(const unsigned int i)
		{
			return m_x[i];
		}

		FORCE_INLINE const Vector3r &getPosition(const unsigned int i) const
		{
			return m_x[i];
		}

		FORCE_INLINE void setPosition(const unsigned int i, const Vector3r &pos)
		{
			fprintf (stderr,"setPosition\n");
			m_x[i] = pos;
		}

		/** Resize the array containing the particle data.
		*/
		FORCE_INLINE void resize(const unsigned int newSize)
		{
			m_x.resize(newSize);
		}

		/** Reserve the array containing the particle data.
		*/
		FORCE_INLINE void reserve(const unsigned int newSize)
		{
			m_x.reserve(newSize);
		}

		/** Release the array containing the particle data.
		*/
		FORCE_INLINE void release()
		{
			m_x.clear();
		}

		/** Release the array containing the particle data.
		*/
		FORCE_INLINE unsigned int size() const
		{
			return (unsigned int)m_x.size();
		}

		FORCE_INLINE const std::vector<Vector3r>* getVertices()
		{
			return &m_x;
		}
	};

	/** This class encapsulates the state of all vertices of a current network model.
	 * All parameters are stored in individual arrays.
	 */
	class CurrentData
	{
		private:
			// 'current network' vertices
			std::vector<Vector3r> m_currA;
			std::vector<Vector3r> m_currB;

		public:
			FORCE_INLINE CurrentData(void)	:
				  m_currA(),
				  m_currB()
			{
			}

			FORCE_INLINE ~CurrentData(void) 
			{
				m_currA.clear();
				m_currB.clear();
			}

			FORCE_INLINE void addVertex(void)
			{
				m_currA.push_back(Vector3r(0.0, 0.0, 0.0));
				m_currB.push_back(Vector3r(0.0, 0.0, 0.0));
			}

			FORCE_INLINE const Vector3r getCurrA(const unsigned int i) const
			{
				return m_currA[i];
			}

			FORCE_INLINE Vector3r& getCurrA(unsigned int i)
			{
				return m_currA[i];
			}

			FORCE_INLINE void setCurrA(unsigned int i, Vector3r inter)
			{
				m_currA[i] = inter;
			}

			FORCE_INLINE const Vector3r getCurrB(const unsigned int i) const
			{
				return m_currB[i];
			}

			FORCE_INLINE Vector3r& getCurrB(unsigned int i)
			{
				return m_currB[i];
			}

			FORCE_INLINE void setCurrB(unsigned int i, Vector3r inter)
			{
				m_currB[i] = inter;
			}

			/** Release the array containing the current network data.
			 */
			FORCE_INLINE void release()
			{
				m_currA.clear();
				m_currB.clear();
			}

			FORCE_INLINE unsigned int size() const 
			{
				return (unsigned int) m_currA.size();
			}
	};

	/** This class encapsulates the state of all particles of a particle model.
	 * All parameters are stored in individual arrays.
	 */
	class ParticleData
	{
		private:
			// Mass
			// If the mass is zero, the particle is static
			std::vector<Real> m_masses;
			std::vector<Real> m_invMasses;

			// radius of particle
			std::vector<Real> m_radius;

			// intersections (up to three pairs: [0,1] [2,3] [4,5])
			// std::vector<Vector3r> m_int[6];
			std::vector<Vector3r> m_int0;
			std::vector<Vector3r> m_int1;
			std::vector<Vector3r> m_int2;
			std::vector<Vector3r> m_int3;
			std::vector<Vector3r> m_int4;
			std::vector<Vector3r> m_int5;
			std::vector<Vector3r> m_int6;
			std::vector<Vector3r> m_int7;
			// Vector2r	m_int;
			std::vector<Real> m_nint;

			// contact with other circles (spheres) - up to three
			std::vector<Real> m_which0;
			std::vector<Real> m_which1;
			std::vector<Real> m_which2;
			std::vector<Real> m_which3;
			std::vector<Real> m_which4;
			std::vector<Real> m_which5;
			std::vector<Real> m_which6;
			std::vector<Real> m_which7;
			std::vector<Real> m_which8;

			// Dynamic state
			std::vector<Vector3r> m_x0;
			std::vector<Vector3r> m_x;
			std::vector<Vector3r> m_v;
			std::vector<Vector3r> m_a;
			std::vector<Vector3r> m_oldX;
			std::vector<Vector3r> m_lastX;

		public:
			FORCE_INLINE ParticleData(void)	:
				  m_masses(),
				  m_invMasses(),
				  m_radius(),
				  m_int0(),
				  m_int1(),
				  m_int2(),
				  m_int3(),
				  m_int4(),
				  m_int5(),
				  m_int6(),
				  m_int7(),
				  m_nint(),
				  m_which0(),
				  m_which1(),
				  m_which2(),
				  m_which3(),
				  m_which4(),
				  m_which5(),
				  m_which6(),
				  m_which7(),
				  m_which8(),
				  m_x0(),
				  m_x(),
				  m_v(),
				  m_a(),
				  m_oldX(),
				  m_lastX()
			{
			}

			FORCE_INLINE ~ParticleData(void) 
			{
				m_masses.clear();
				m_invMasses.clear();
				m_radius.clear();
				m_int0.clear();
				m_int1.clear();
				m_int2.clear();
				m_int3.clear();
				m_int4.clear();
				m_int5.clear();
				m_int6.clear();
				m_int7.clear();
				m_nint.clear();
				m_which0.clear();
				m_which1.clear();
				m_which2.clear();
				m_which3.clear();
				m_which4.clear();
				m_which5.clear();
				m_which6.clear();
				m_which7.clear();
				m_which8.clear();

				m_x0.clear();
				m_x.clear();
				m_v.clear();
				m_a.clear();
				m_oldX.clear();
				m_lastX.clear();
			}

			FORCE_INLINE void addVertex(const Vector3r &vertex)
			{
				m_x0.push_back(vertex);
				m_x.push_back(vertex);
				m_oldX.push_back(vertex);
				m_lastX.push_back(vertex);
				m_masses.push_back(1.0);
				m_invMasses.push_back(1.0);

				m_radius.push_back(0.0);

				m_int0.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int1.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int2.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int3.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int4.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int5.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int6.push_back(Vector3r(0.0, 0.0, 0.0));
				m_int7.push_back(Vector3r(0.0, 0.0, 0.0));
				m_nint.push_back(0);

				m_which0.push_back(0);
				m_which1.push_back(0);
				m_which2.push_back(0);
				m_which3.push_back(0);
				m_which4.push_back(0);
				m_which5.push_back(0);
				m_which6.push_back(0);
				m_which7.push_back(0);
				m_which8.push_back(0);

				m_v.push_back(Vector3r(0.0, 0.0, 0.0));
				m_a.push_back(Vector3r(0.0, 0.0, 0.0));
			}

			FORCE_INLINE Vector3r &getPosition(const unsigned int i)
			{
				return m_x[i];
			}

			FORCE_INLINE const Vector3r &getPosition(const unsigned int i) const 
			{
				return m_x[i];
			}

			FORCE_INLINE void setPosition(const unsigned int i, const Vector3r &pos)
			{
				m_x[i] = pos;
			}

			FORCE_INLINE Vector3r &getPosition0(const unsigned int i)
			{
				return m_x0[i];
			}

			FORCE_INLINE const Vector3r &getPosition0(const unsigned int i) const
			{
				return m_x0[i];
			}

			FORCE_INLINE void setPosition0(const unsigned int i, const Vector3r &pos)
			{
				m_x0[i] = pos;
			}

			FORCE_INLINE Vector3r &getLastPosition(const unsigned int i)
			{
				return m_lastX[i];
			}

			FORCE_INLINE const Vector3r &getLastPosition(const unsigned int i) const
			{
				return m_lastX[i];
			}

			FORCE_INLINE void setLastPosition(const unsigned int i, const Vector3r &pos)
			{
				m_lastX[i] = pos;
			}

			FORCE_INLINE Vector3r &getOldPosition(const unsigned int i)
			{
				return m_oldX[i];
			}

			FORCE_INLINE const Vector3r &getOldPosition(const unsigned int i) const
			{
				return m_oldX[i];
			}

			FORCE_INLINE void setOldPosition(const unsigned int i, const Vector3r &pos)
			{
				m_oldX[i] = pos;
			}
			
			FORCE_INLINE Vector3r &getVelocity(const unsigned int i)
			{
				return m_v[i];
			}

			FORCE_INLINE const Vector3r &getVelocity(const unsigned int i) const 
			{
				return m_v[i];
			}

			FORCE_INLINE void setVelocity(const unsigned int i, const Vector3r &vel)
			{
				m_v[i] = vel;
			}

			FORCE_INLINE Vector3r &getAcceleration(const unsigned int i)
			{
				return m_a[i];
			}

			FORCE_INLINE const Vector3r &getAcceleration(const unsigned int i) const 
			{
				return m_a[i];
			}

			FORCE_INLINE void setAcceleration(const unsigned int i, const Vector3r &accel)
			{
				m_a[i] = accel;
			}

			FORCE_INLINE const Real getMass(const unsigned int i) const
			{
				return m_masses[i];
			}

			FORCE_INLINE Real& getMass(const unsigned int i)
			{
				return m_masses[i];
			}

			FORCE_INLINE void setMass(const unsigned int i, const Real mass)
			{
				m_masses[i] = mass;
				if (mass != 0.0)
					m_invMasses[i] = 1.0 / mass;
				else
					m_invMasses[i] = 0.0;
			}

			FORCE_INLINE const Real getRadius(const unsigned int i) const
			{
				return m_radius[i];
			}

			FORCE_INLINE Real& getRadius(const unsigned int i)
			{
				return m_radius[i];
			}

			FORCE_INLINE void setRadius(unsigned int i, Real radius)
			{
				m_radius[i] = radius;
			}

			FORCE_INLINE const Real getWhich(const unsigned int i, const unsigned int which) const
			{
				if (which == 0) {
					return m_which0[i];
				} else if (which == 1) {
					return m_which1[i];
				} else if (which == 2) {
					return m_which2[i];
				} else if (which == 3) {
					return m_which3[i];
				} else if (which == 4) {
					return m_which4[i];
				} else if (which == 5) {
					return m_which5[i];
				} else if (which == 6) {
					return m_which6[i];
				} else if (which == 7) {
					return m_which7[i];
				} else if (which == 8) {
					return m_which8[i];
				}
			}

			FORCE_INLINE Real& getWhich(unsigned int i, unsigned int which)
			{
				if (which == 0) {
					return m_which0[i];
				} else if (which == 1) {
					return m_which1[i];
				} else if (which == 2) {
					return m_which2[i];
				} else if (which == 3) {
					return m_which3[i];
				} else if (which == 4) {
					return m_which4[i];
				} else if (which == 5) {
					return m_which5[i];
				} else if (which == 6) {
					return m_which6[i];
				} else if (which == 7) {
					return m_which7[i];
				} else if (which == 8) {
					return m_which8[i];
				}
			}

			FORCE_INLINE const Vector3r getInt(const unsigned int i, const unsigned int which) const
			{
				if (which == 0) {
					return m_int0[i];
				} else if (which == 1) {
					return m_int1[i];
				} else if (which == 2) {
					return m_int2[i];
				} else if (which == 3) {
					return m_int3[i];
				} else if (which == 4) {
					return m_int4[i];
				} else if (which == 5) {
					return m_int5[i];
				} else if (which == 6) {
					return m_int6[i];
				} else if (which == 7) {
					return m_int7[i];
				}
				// return m_int;
			}

			FORCE_INLINE Vector3r& getInt(unsigned int i, unsigned int which)
			{
				if (which == 0) {
					return m_int0[i];
				} else if (which == 1) {
					return m_int1[i];
				} else if (which == 2) {
					return m_int2[i];
				} else if (which == 3) {
					return m_int3[i];
				} else if (which == 4) {
					return m_int4[i];
				} else if (which == 5) {
					return m_int5[i];
				} else if (which == 6) {
					return m_int6[i];
				} else if (which == 7) {
					return m_int7[i];
				}
				// return m_int;
			}

			FORCE_INLINE void setInt(unsigned int i, Vector3r inter, unsigned int which)
			{
				if (which == 0) {
					m_int0[i] = inter;
				} else if (which == 1) {
					m_int1[i] = inter;
				} else if (which == 2) {
					m_int2[i] = inter;
				} else if (which == 3) {
					m_int3[i] = inter;
				} else if (which == 4) {
					m_int4[i] = inter;
				} else if (which == 5) {
					m_int5[i] = inter;
				} else if (which == 6) {
					m_int6[i] = inter;
				} else if (which == 7) {
					m_int7[i] = inter;
				}
				// m_int = inter;
			}

			FORCE_INLINE const Real getnInt(const unsigned int i) const
			{
				return m_nint[i];
			}

			FORCE_INLINE Real& getnInt(const unsigned int i)
			{
				return m_nint[i];
			}

			FORCE_INLINE void addInt(unsigned int i, unsigned int which, unsigned int n)
			{
				fprintf (stderr,"addInt %d %d %d\n",i,which,n);
				// sets [which] circle is in contact with current [i] circle
				if (n == 0) {
					m_which0[i] = which;
				} else if (n == 1) {
					m_which1[i] = which;
				} else if (n == 2) {
					m_which2[i] = which;
				} else if (n == 3) {
					m_which3[i] = which;
				} else if (n == 4) {
					m_which4[i] = which;
				} else if (n == 5) {
					m_which5[i] = which;
				} else if (n == 6) {
					m_which6[i] = which;
				} else if (n == 7) {
					m_which7[i] = which;
				} else if (n == 8) {
					m_which8[i] = which;
				}
				m_nint[i] += (Real)1;
			}

			FORCE_INLINE void setnInt(unsigned int i, Real inter)
			{
				m_nint[i] = inter;
			}

			FORCE_INLINE const Real getInvMass(const unsigned int i) const
			{
				return m_invMasses[i];
			}

			FORCE_INLINE const unsigned int getNumberOfParticles() const
			{
				return (unsigned int) m_x.size();
			}

			/** Resize the array containing the particle data.
			 */
			FORCE_INLINE void resize(const unsigned int newSize)
			{
				m_masses.resize(newSize);
				m_invMasses.resize(newSize);
				m_radius.resize(newSize);
				m_x0.resize(newSize);
				m_x.resize(newSize);
				m_v.resize(newSize);
				m_a.resize(newSize);
				m_oldX.resize(newSize);
				m_lastX.resize(newSize);
			}

			/** Reserve the array containing the particle data.
			 */
			FORCE_INLINE void reserve(const unsigned int newSize)
			{
				m_masses.reserve(newSize);
				m_invMasses.reserve(newSize);
				m_radius.reserve(newSize);
				m_x0.reserve(newSize);
				m_x.reserve(newSize);
				m_v.reserve(newSize);
				m_a.reserve(newSize);
				m_oldX.reserve(newSize);
				m_lastX.reserve(newSize);
			}

			/** Release the array containing the particle data.
			 */
			FORCE_INLINE void release()
			{
				m_masses.clear();
				m_invMasses.clear();
				m_radius.clear();
				m_x0.clear();
				m_x.clear();
				m_v.clear();
				m_a.clear();
				m_oldX.clear();
				m_lastX.clear();
			}

			/** Release the array containing the particle data.
			 */
			FORCE_INLINE unsigned int size() const 
			{
				return (unsigned int) m_x.size();
			}
	};

	/** This class encapsulates the state of all orientations of a quaternion model.
	* All parameters are stored in individual arrays.
	*/
	class OrientationData
	{
	private:
		// Mass
		// If the mass is zero, the particle is static
		std::vector<Real> m_masses;
		std::vector<Real> m_invMasses;

		// Dynamic state
		std::vector<Quaternionr> m_q0;
		std::vector<Quaternionr> m_q;
		std::vector<Vector3r> m_omega;
		std::vector<Vector3r> m_alpha;
		std::vector<Quaternionr> m_oldQ;
		std::vector<Quaternionr> m_lastQ;

	public:
		FORCE_INLINE OrientationData(void) :
			m_masses(),
			m_invMasses(),
			m_q0(),
			m_q(),
			m_omega(),
			m_alpha(),
			m_oldQ(),
			m_lastQ()
		{
		}

		FORCE_INLINE ~OrientationData(void)
		{
			m_masses.clear();
			m_invMasses.clear();
			m_q0.clear();
			m_q.clear();
			m_omega.clear();
			m_alpha.clear();
			m_oldQ.clear();
			m_lastQ.clear();
		}

		FORCE_INLINE void addQuaternion(const Quaternionr &vertex)
		{
			m_q0.push_back(vertex);
			m_q.push_back(vertex);
			m_oldQ.push_back(vertex);
			m_lastQ.push_back(vertex);
			m_masses.push_back(1.0);
			m_invMasses.push_back(1.0);
			m_omega.push_back(Vector3r(0.0, 0.0, 0.0));
			m_alpha.push_back(Vector3r(0.0, 0.0, 0.0));
		}

		FORCE_INLINE Quaternionr &getQuaternion(const unsigned int i)
		{
			return m_q[i];
		}

		FORCE_INLINE const Quaternionr &getQuaternion(const unsigned int i) const
		{
			return m_q[i];
		}

		FORCE_INLINE void setQuaternion(const unsigned int i, const Quaternionr &pos)
		{
			m_q[i] = pos;
		}

		FORCE_INLINE Quaternionr &getQuaternion0(const unsigned int i)
		{
			return m_q0[i];
		}

		FORCE_INLINE const Quaternionr &getQuaternion0(const unsigned int i) const
		{
			return m_q0[i];
		}

		FORCE_INLINE void setQuaternion0(const unsigned int i, const Quaternionr &pos)
		{
			m_q0[i] = pos;
		}

		FORCE_INLINE Quaternionr &getLastQuaternion(const unsigned int i)
		{
			return m_lastQ[i];
		}

		FORCE_INLINE const Quaternionr &getLastQuaternion(const unsigned int i) const
		{
			return m_lastQ[i];
		}

		FORCE_INLINE void setLastQuaternion(const unsigned int i, const Quaternionr &pos)
		{
			m_lastQ[i] = pos;
		}

		FORCE_INLINE Quaternionr &getOldQuaternion(const unsigned int i)
		{
			return m_oldQ[i];
		}

		FORCE_INLINE const Quaternionr &getOldQuaternion(const unsigned int i) const
		{
			return m_oldQ[i];
		}

		FORCE_INLINE void setOldQuaternion(const unsigned int i, const Quaternionr &pos)
		{
			m_oldQ[i] = pos;
		}

		FORCE_INLINE Vector3r &getVelocity(const unsigned int i)
		{
			return m_omega[i];
		}

		FORCE_INLINE const Vector3r &getVelocity(const unsigned int i) const
		{
			return m_omega[i];
		}

		FORCE_INLINE void setVelocity(const unsigned int i, const Vector3r &vel)
		{
			m_omega[i] = vel;
		}

		FORCE_INLINE Vector3r &getAcceleration(const unsigned int i)
		{
			return m_alpha[i];
		}

		FORCE_INLINE const Vector3r &getAcceleration(const unsigned int i) const
		{
			return m_alpha[i];
		}

		FORCE_INLINE void setAcceleration(const unsigned int i, const Vector3r &accel)
		{
			m_alpha[i] = accel;
		}

		FORCE_INLINE const Real getMass(const unsigned int i) const
		{
			return m_masses[i];
		}

		FORCE_INLINE Real& getMass(const unsigned int i)
		{
			return m_masses[i];
		}

		FORCE_INLINE void setMass(const unsigned int i, const Real mass)
		{
			m_masses[i] = mass;
			if (mass != 0.0)
				m_invMasses[i] = 1.0 / mass;
			else
				m_invMasses[i] = 0.0;
		}

		FORCE_INLINE const Real getInvMass(const unsigned int i) const
		{
			return m_invMasses[i];
		}

		FORCE_INLINE const unsigned int getNumberOfQuaternions() const
		{
			return (unsigned int)m_q.size();
		}

		/** Resize the array containing the particle data.
		*/
		FORCE_INLINE void resize(const unsigned int newSize)
		{
			m_masses.resize(newSize);
			m_invMasses.resize(newSize);
			m_q0.resize(newSize);
			m_q.resize(newSize);
			m_omega.resize(newSize);
			m_alpha.resize(newSize);
			m_oldQ.resize(newSize);
			m_lastQ.resize(newSize);
		}

		/** Reserve the array containing the particle data.
		*/
		FORCE_INLINE void reserve(const unsigned int newSize)
		{
			m_masses.reserve(newSize);
			m_invMasses.reserve(newSize);
			m_q0.reserve(newSize);
			m_q.reserve(newSize);
			m_omega.reserve(newSize);
			m_alpha.reserve(newSize);
			m_oldQ.reserve(newSize);
			m_lastQ.reserve(newSize);
		}

		/** Release the array containing the particle data.
		*/
		FORCE_INLINE void release()
		{
			m_masses.clear();
			m_invMasses.clear();
			m_q0.clear();
			m_q.clear();
			m_omega.clear();
			m_alpha.clear();
			m_oldQ.clear();
			m_lastQ.clear();
		}

		/** Release the array containing the particle data.
		*/
		FORCE_INLINE unsigned int size() const
		{
			return (unsigned int)m_q.size();
		}
	};
}

#endif
