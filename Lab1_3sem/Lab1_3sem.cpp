#include <iostream>
#include <cmath>
#include <string> 
#include <numbers>

#define CIR_GEOMETRY 3
#define RECT_GEOMETRY 4
#define R 5
#define SIDE_MYLTIPLIER 10

using namespace std;

double get_PI() {
    return std::numbers::pi;
}

bool cir_or_rect(double probability = 0.5) {
    return (rand() / (double)(RAND_MAX)) < probability;
}

double angle_rand() {
    return (double)(rand()) / RAND_MAX * (get_PI() / 2);
}

double rand_in_range(double min, double max) {
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}

class Figure {
public:
    static int numberOfFigures;
    double* geometry;
    double angle;

    // Чисто виртуальная функция (абстрактная функция)
    virtual void rotate() = 0;
    virtual double square() = 0;
    virtual string toString() = 0;
    // Виртуальный деструктор
    virtual ~Figure() = default;
};

int Figure::numberOfFigures = 0;

class Circle : public Figure {
public:
    double* geometry = new double[CIR_GEOMETRY]; // xp yp r

    Circle() {
        geometry[0] = 0.0; 
        geometry[1] = 0.0; 
        geometry[2] = 1.0; 
        angle = 0.0; 
        Figure::numberOfFigures++;
    }
    Circle(double xp, double yp, double r, double angle) {
        geometry[0] = xp, geometry[1] = yp, geometry[2] = r;
        this->angle = angle;
        Figure::numberOfFigures++;
    }

    void rotate() override {}

    double square() override {
        return get_PI() * geometry[2] * geometry[2];
    }

    string toString() override {
        return string(typeid(*this).name()) +
            " radius:" + to_string(geometry[2]) +
            " c_point(" + to_string(geometry[0]) + ";" + to_string(geometry[1]) + ")";
    }

    static bool checkCircleInBox(Circle circle) {

        double x = circle.geometry[0];
        double y = circle.geometry[1];
        double r = circle.geometry[2]; // && y_point > r && 10*R - x_point > r && 10*R - y_point > r

        return (x - r >= 0 && y - r >= 0 && x + r <= 10 * R && y + r <= 10 * R);
    }

    bool checkIsInBox() {

        double x = geometry[0];
        double y = geometry[1];
        double r = geometry[2]; // && y_point > r && 10*R - x_point > r && 10*R - y_point > r

        return (x - r >= 0 && y - r >= 0 && x + r <= 10 * R && y + r <= 10 * R);
    }

    ~Circle() {
        delete[] geometry;
        Figure::numberOfFigures--;
    }
};

class Rectangle : public Figure {
public:
    double* geometry = new double[RECT_GEOMETRY]; // xp yp a b

    Rectangle() {
        geometry[0] = 0.0; // x координата центра
        geometry[1] = 0.0; // y координата центра
        geometry[2] = 1.0; // ширина (a), например 1
        geometry[3] = 1.0; // высота (b), например 1
        angle = 0.0; // угол
        Figure::numberOfFigures++;
    }
    Rectangle(double xp, double yp, double a, double b, double angle) {
        geometry[0] = xp, geometry[1] = yp, geometry[2] = a, geometry[3] = b;
        this->angle = angle;
        Figure::numberOfFigures++;
    }

    void rotate() override {}

    double square() override {
        return geometry[2] * geometry[3];
    }

    string toString() override {
        return string(typeid(*this).name()) +
            " sides(a;b):(" + to_string(geometry[2]) + ";" + to_string(geometry[3]) + ")" +
            " c_point(" + to_string(geometry[0]) + ";" + to_string(geometry[1]) + ")";
    }

    static bool checkRectangleInBox(Rectangle rectangle) {

    double x = rectangle.geometry[0]; // x координата центра
    double y = rectangle.geometry[1]; // y координата центра
    double a = rectangle.geometry[2]; // ширина прямоугольника
    double b = rectangle.geometry[3]; // высота прямоугольника
    double angleRad = rectangle.angle * get_PI() / 180.0; // угол в радианах


    double half_a = a / 2;
    double half_b = b / 2;

    double corners[4][2] = {
        {-half_a, -half_b}, // нижний левый
        { half_a, -half_b}, // нижний правый
        { half_a,  half_b}, // верхний правый
        {-half_a,  half_b}  // верхний левый
    };

    for (int i = 0; i < 4; ++i) {
        double local_x = corners[i][0];
        double local_y = corners[i][1];

        // Применяем матрицу поворота
        double rotated_x = x + (local_x * cos(angleRad) - local_y * sin(angleRad));
        double rotated_y = y + (local_x * sin(angleRad) + local_y * cos(angleRad));

        // Проверяем, что угол находится внутри границ
        if (rotated_x < 0 || rotated_x > 30 || rotated_y < 0 || rotated_y > 30) {
            return false; // Если хотя бы один угол выходит за пределы, возвращаем false
        }
    }

    return true; // Все углы внутри коробки
}

    bool checkIsInBox() {
        double x = geometry[0]; // x координата центра
        double y = geometry[1]; // y координата центра
        double a = geometry[2]; // ширина прямоугольника
        double b = geometry[3]; // высота прямоугольника
        double angleRad = angle * get_PI() / 180.0; // угол в радианах


        double half_a = a / 2;
        double half_b = b / 2;

        double corners[4][2] = {
            {-half_a, -half_b}, // нижний левый
            { half_a, -half_b}, // нижний правый
            { half_a,  half_b}, // верхний правый
            {-half_a,  half_b}  // верхний левый
        };

        for (int i = 0; i < 4; ++i) {
            double local_x = corners[i][0];
            double local_y = corners[i][1];

            // Применяем матрицу поворота
            double rotated_x = x + (local_x * cos(angleRad) - local_y * sin(angleRad));
            double rotated_y = y + (local_x * sin(angleRad) + local_y * cos(angleRad));

            // Проверяем, что угол находится внутри границ
            if (rotated_x < 0 || rotated_x > 30 || rotated_y < 0 || rotated_y > 30) {
                return false; // Если хотя бы один угол выходит за пределы, возвращаем false
            }
        }

        return true; // Все углы внутри коробки
    }

    ~Rectangle() {
        delete[] geometry;
        Figure::numberOfFigures--;
    }
};

std::ostream& operator << (std::ostream& stream, Figure& c1) {
    stream << "Value of class: ";
    stream << c1.toString();
    return stream;
}

class Field {
private:

    void generate_common_params(double& xp, double& yp, double& angle) {
        xp = rand_in_range(0, 10 * R);
        yp = rand_in_range(0, 10 * R);
        angle = angle_rand();
    }

    bool create_circle(Figure** exeplars, int i) {
        double xp, yp, r, angle;
        generate_common_params(xp, yp, angle);
        r = rand_in_range(0.1 * R, 0.5 * R);

        Circle circle(xp, yp, r, angle);
        if (!circle.checkIsInBox()) {
            return false; // не попал в коробку
        }

        exeplars[i] = new Circle(xp, yp, r, angle); // создаем объект
        return true;
    }

    bool create_rectangle(Figure** exeplars, int i) {
        double xp, yp, a, b, angle;
        generate_common_params(xp, yp, angle);
        a = rand_in_range(0.1 * R, 0.5 * R);
        b = rand_in_range(0.1 * R, 0.5 * R);

        Rectangle rectangle(xp, yp, a, b, angle);
        if (!rectangle.checkIsInBox()) {
            return false; 
        }

        exeplars[i] = new Rectangle(xp, yp, a, b, angle);
        return true;
    }

    Figure** generate_figures(const int K = 1) {
        Figure** exeplars = new Figure * [K];
        for (int i = 0; i < K; i++) {
            if (cir_or_rect(probabilityCircle)) {
                if (!create_circle(exeplars, i)) {
                    i--; // если не попал 
                }
            }
            else {
                if (!create_rectangle(exeplars, i)) {
                    i--; // если не попал 
                }
            }
        }
        return exeplars;
    }

    bool circlesIntersect(Circle* c1, Circle* c2) {
        double dx = c1->geometry[0] - c2->geometry[0]; // Разница по X
        double dy = c1->geometry[1] - c2->geometry[1]; // Разница по Y
        double distance = sqrt(dx * dx + dy * dy); // Расстояние между центрами
        double rSum = c1->geometry[2] + c2->geometry[2]; // Сумма радиусов

        return distance <= rSum; // Если расстояние меньше суммы радиусов, круги пересекаются
    }

    bool rectangleCircleIntersect(Rectangle* rect, Circle* circle) {
        double cx = circle->geometry[0];
        double cy = circle->geometry[1];
        double r = circle->geometry[2];

        double rectX = rect->geometry[0];
        double rectY = rect->geometry[1];
        double a = rect->geometry[2] / 2; // Половина ширины
        double b = rect->geometry[3] / 2; // Половина высоты

        double closestX = std::max(rectX - a, std::min(cx, rectX + a));
        double closestY = std::max(rectY - b, std::min(cy, rectY + b));

        double dx = cx - closestX;
        double dy = cy - closestY;

        return (dx * dx + dy * dy) <= r * r; // Если расстояние меньше радиуса, то они пересекаются
    }

    bool rectanglesIntersect(Rectangle* r1, Rectangle* r2) {
        double r1_left = r1->geometry[0] - r1->geometry[2] / 2;
        double r1_right = r1->geometry[0] + r1->geometry[2] / 2;
        double r1_top = r1->geometry[1] + r1->geometry[3] / 2;
        double r1_bottom = r1->geometry[1] - r1->geometry[3] / 2;

        double r2_left = r2->geometry[0] - r2->geometry[2] / 2;
        double r2_right = r2->geometry[0] + r2->geometry[2] / 2;
        double r2_top = r2->geometry[1] + r2->geometry[3] / 2;
        double r2_bottom = r2->geometry[1] - r2->geometry[3] / 2;

        // Прямоугольники пересекаются, если одна из сторон одного прямоугольника пересекает другой
        return !(r1_left > r2_right || r1_right < r2_left || r1_top < r2_bottom || r1_bottom > r2_top);
    }

    bool figuresIntersect(Figure* f1, Figure* f2) {
        if (typeid(*f1) == typeid(Circle) && typeid(*f2) == typeid(Circle)) {
            // Круг-круг
            return circlesIntersect(static_cast<Circle*>(f1), static_cast<Circle*>(f2));
        }
        else if (typeid(*f1) == typeid(Rectangle) && typeid(*f2) == typeid(Rectangle)) {
            // Прямоугольник-прямоугольник
            return rectanglesIntersect(static_cast<Rectangle*>(f1), static_cast<Rectangle*>(f2));
        }
        else if (typeid(*f1) == typeid(Circle) && typeid(*f2) == typeid(Rectangle)) {
            // Круг-прямоугольник
            return rectangleCircleIntersect(static_cast<Rectangle*>(f2), static_cast<Circle*>(f1));
        }
        else if (typeid(*f1) == typeid(Rectangle) && typeid(*f2) == typeid(Circle)) {
            // Прямоугольник-круг
            return rectangleCircleIntersect(static_cast<Rectangle*>(f1), static_cast<Circle*>(f2));
        }
        return false; // Если типы не поддерживаются
    }

    Figure* create_combined_figure(Figure* fig1, Figure* fig2) {
        Figure* larger = (fig1->square() > fig2->square()) ? fig1 : fig2;
        double new_area = fig1->square() + fig2->square();

        if (Circle* larger_circle = dynamic_cast<Circle*>(larger)) {
            double new_radius = sqrt(new_area / get_PI()); // Площадь круга = πr² => r = sqrt(area / π)
            return new Circle(larger_circle->geometry[0], larger_circle->geometry[1], new_radius, larger_circle->angle);
        }
        else if (Rectangle* larger_rect = dynamic_cast<Rectangle*>(larger)) {
            double aspect_ratio = larger_rect->geometry[2] / larger_rect->geometry[3]; // Сохранение пропорций
            double new_a = sqrt(new_area * aspect_ratio);
            double new_b = new_a / aspect_ratio;
            return new Rectangle(larger_rect->geometry[0], larger_rect->geometry[1], new_a, new_b, larger_rect->angle);
        }
        return nullptr; // Если непонятная фигурка
    }

public:
    Figure** exeplars;
    const double BOX_SIDE = SIDE_MYLTIPLIER * R;
    int figures_amount;
    double probabilityCircle = 0.5;

    void setProbCircle(double probability) {
        probabilityCircle = probability;
    }

    void newField(const int K = 1) {
       exeplars = generate_figures(K);
    }


    void printFigures() {
        for (int i = 0; i < figures_amount; i++ ) {
            cout << *exeplars[i] << " number :" << i << endl;
        }
        cout << "All figures amount: " << Figure::numberOfFigures << endl;;
    }

    void combine(){
        for (int i = 0; i < figures_amount; i++) {
            for (int j = i + 1; j < figures_amount; j++) {
                if (figuresIntersect(exeplars[i], exeplars[j])) {
                    Figure* combined_figure = create_combined_figure(exeplars[i], exeplars[j]);
                    delete exeplars[i];
                    delete exeplars[j];
                    exeplars[i] = combined_figure;
                    for (int k = j; k < figures_amount - 1; k++) {
                        exeplars[k] = exeplars[k + 1];
                    }
                    figures_amount--;
                    j--;
                }
            }
        }
    }
    void squreAllFiguresPrint() {
        double resultSquare = 0;
        for (int i = 0; i < figures_amount; i++) {
            resultSquare += exeplars[i]->square();
        }
        cout << "All figures squares: " << resultSquare << endl;
    }

    Field(int amount_of_figures, double probability_of_circle = 0.5) {
        figures_amount = amount_of_figures;
        setProbCircle(probability_of_circle);
        exeplars = generate_figures(amount_of_figures);
        
    }

    ~Field() {
        for (int i = 0; i < figures_amount; i++) {
            delete exeplars[i];  
        }
        delete[] exeplars;  
    }

};




int main() {
    Field field = Field(20);
    field.printFigures();
    field.squreAllFiguresPrint();
    field.combine();
    field.printFigures();
    field.squreAllFiguresPrint();
}

