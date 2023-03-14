#include "matrix4x4.h"
#include "quaternion.h"

void Matrix4x4::rotate(const Vector3& rotate) {
	Quaternion q = Quaternion::from_euler(rotate);
	Matrix4x4 m = q.to_matrix4x4();
	*this *= m;
}
