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
    virtual void rotate() = 0;
    virtual void square() = 0;
    virtual string toString() = 0;
    // Виртуальный деструктор
    virtual ~Figure() = default;
};


class Circle : public Figure {
public:
    double* geometry = new double[CIR_GEOMETRY]; // xp yp r
    Circle();
    Circle(double xp, double yp, double r, double angle) {
        geometry[0] = xp, geometry[1] = yp, geometry[2] = r;
        this->angle = angle;
    }
    void rotate() override {}
    void square() override {}
    string toString() override {
        string result = typeid(this).name()
            + string(" radius:") + to_string(geometry[2])
            + string(" c_point(") + to_string(geometry[0]) + string(";") + to_string(geometry[1]) + string(")");
        return result;
    }
    ~Circle() override {
        delete[] geometry;
    }
};


class Rectangle : public Figure {
public:
    double* geometry = new double[RECT_GEOMETRY]; // xp yp a b
    Rectangle();
    Rectangle(double xp, double yp, double a, double b, double angle) {
        geometry[0] = xp, geometry[1] = yp, geometry[2] = a, geometry[3] = b;
        this->angle = angle;
    }
    void rotate() override {}
    void square() override {}
    string toString() override {
        string result = typeid(this).name()
            + string(" sides(a;b):(") + to_string(geometry[2]) + string(";") + to_string(geometry[3]) + string(")")
            + string(" c_point(") + to_string(geometry[0]) + string(";") + to_string(geometry[1]) + string(")");
        return result;
    }
    ~Rectangle() {
        delete[] geometry;
    }
};
double get_PI() {
    double pi = 3.14;
    return pi;
}

bool cir_or_rect() {
    return rand() % 2 == 0 ? true : false;
}

double angle_rand() {
    return (double)(rand()) / RAND_MAX * (get_PI() / 2);
}

double rand_in_range(double min, double max) {
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}
bool kill_circle_or_not(Circle circle) {
    double x_point = circle.geometry[0];
    double y_point = circle.geometry[1];
    double r = circle.geometry[2]; // && y_point > r && 10*R - x_point > r && 10*R - y_point > r
    if (x_point >= r && y_point >= r && 10 * R - x_point > r && 10 * R - y_point > r) {
        return true;
    }
    return false;
}

Figure** generate_figures(const int K = 1) {
    Figure** exeplars = new Figure * [MAS_SIZE];

    for (int i = 0; i < K; i++) {
        if (cir_or_rect()) {
            exeplars[i] = new Circle(
                rand_in_range(0, 10 * R),
                rand_in_range(0, 10 * R),
                rand_in_range(0.1 * R, 0.5 * R),
                angle_rand()
            );
        }
        else {
            exeplars[i] = new Rectangle(
                rand_in_range(0, 10 * R),
                rand_in_range(0, 10 * R),
                rand_in_range(0.1 * R, 0.5 * R),
                rand_in_range(0.1 * R, 0.5 * R),
                angle_rand()
            );
        }
    }
    return exeplars;
}


int main() {
    Figure** exeplars = generate_figures(20);
    for (int i = 0; i < 20; i++) {
        cout << exeplars[i]->toString() << "||" << endl;
    }
}