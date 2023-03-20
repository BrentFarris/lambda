#include <iostream>
#include <identities/entity.h>
#include <identities/transformations.h>

using namespace lambda;

int main() {
	Transformations transformations(10);
	Entity a(transformations);
	Entity b(transformations);
	b.transform.set_parent(a.transform);
	b.transform.set_position(Vector3(0.0F, 1.0F, 0.0F));
	a.transform.set_rotation(Vector3(0.0F, 0.0F, 90.0F));
	transformations.update();

	std::cout << b.transform.matrix().to_string() << std::endl;
	std::cout << "Hello, World!" << std::endl;
	return 0;
}
