#include <iostream>
#include <cmath>
#include <string> 

#define CIR_GEOMETRY 3
#define RECT_GEOMETRY 4
#define R 5
#define MAS_SIZE 100

using namespace std;

class Figure {
public:
    double* geometry;
    double angle;

    // Чисто виртуальная функция (абстрактная функция)
    virtual void rotate()  = 0;
    virtual void square()  = 0; 
    // Виртуальный деструктор
    virtual ~Figure() = default;
};


class Circle : public Figure {
public:
    double* geometry = new double[CIR_GEOMETRY]; // xp yp r
    Circle();
    Circle(double xp, double yp, double r){
        geometry[0] = xp, geometry[1] = yp, geometry[2] = r;
        this->angle = angle;
    }
    void rotate() override{}
    void square() override{}
};


class Rectangle : public Figure {
public:
    double* geometry = new double[RECT_GEOMETRY]; // xp yp a b
    Rectangle() = default;
    Rectangle(double xp, double yp, double a, double b){
        geometry[0] = xp, geometry[1] = yp, geometry[2] = a, geometry[3] = b;
    }
    void rotate() override{}
    void square() override{}
};


int main(){
    Circle cicrc = Circle(1.0,1.0,1.0);
    Rectangle rect = Rectangle(1.0,1.0,1.0,1.0);
    cicrc.angle = 121;
    rect.angle = 33;
    cout << cicrc.angle << "|" << rect.angle << endl;
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
