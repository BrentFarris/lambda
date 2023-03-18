#include <iostream>
#include <identities/transformations.h>

int main() {
	Transformations transformations(10);
	TransformClaim child(transformations);
	TransformClaim parent(transformations);
	child.set_parent(parent);
	child.set_position(Vector3(0.0F, 1.0F, 0.0F));
	parent.set_rotation(Vector3(0.0F, 0.0F, 90.0F));
	transformations.update();

	std::cout << child.matrix().to_string() << std::endl;
	std::cout << "Hello, World!" << std::endl;
	return 0;
}
