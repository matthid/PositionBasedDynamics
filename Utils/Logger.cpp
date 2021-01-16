#include "Logger.h"

std::ostream& operator<<(std::ostream& os, const Quaternionr& q)
{
	os << "Quaternion(x: " << q.x() << ", y: " << q.y() << ", z: " << q.z() << ", w: " << q.w() << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, const Vector3r& v)
{
	os << "Vector(x: " << v.x() << ", y: " << v.y() << ", z: " << v.z() << ")";
	return os;
}
