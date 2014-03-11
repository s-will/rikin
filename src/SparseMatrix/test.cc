#include <iostream>

class A {
public:
    void
    printA() const {
	std::cout << "A"<< std::endl;
    }
};


class B {
public:
    operator A() {
	return A();
    }

    void
    printB() const {
	std::cout << "B"<< std::endl;
    }
};

void 
print(const A &a){
    a.printA();
}

int 
main() {

    class B b;

    //b.printA();

    print(b);

}
