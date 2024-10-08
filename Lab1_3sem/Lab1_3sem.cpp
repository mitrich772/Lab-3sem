#include <iostream>
#include <cmath>
#include <string> 
#include <numbers>
#include <array>

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
    static inline int numberOfFigures = 0;
    double* geometry;
    double angle;

    // Чисто виртуальная функция (абстрактная функция)
    virtual void rotate(double angle) = 0;
    virtual double square() = 0;
    virtual string toString() = 0;
    // Виртуальный деструктор
    virtual ~Figure() = default;
};

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

    void rotate(double angle) override {
        this->angle = angle;
    }

    double square() override {
        return get_PI() * geometry[2] * geometry[2];
    }

    string toString() override {
        return string(typeid(*this).name()) +
            " radius:" + to_string(geometry[2]) +
            " c_point(" + to_string(geometry[0]) + ";" + to_string(geometry[1]) + ")";
    }

    static bool checkCircleInBox(Circle circle) {
        return circle.checkIsInBox();
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

    void rotate(double angle) override {
        this->angle = angle;
    }

    double square() override {
        return geometry[2] * geometry[3];
    }

    string toString() override {
        return string(typeid(*this).name()) +
            " sides(a;b):(" + to_string(geometry[2]) + ";" + to_string(geometry[3]) + ")" +
            " c_point(" + to_string(geometry[0]) + ";" + to_string(geometry[1]) + ")";
    }

    static bool checkRectangleInBox(Rectangle rectangle) {

        return rectangle.checkIsInBox();
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

    void getCorners(Rectangle* rect, double corners[4][2]) {
        double x = rect->geometry[0]; // центр
        double y = rect->geometry[1]; // центр
        double a = rect->geometry[2] / 2; // половина ширины
        double b = rect->geometry[3] / 2; // половина высоты
        double angleRad = rect->angle * get_PI() / 180.0; // угол в радианах

        // Вершины без учета поворота
        double localCorners[4][2] = {
            {-a, -b}, // нижний левый
            { a, -b}, // нижний правый
            { a,  b}, // верхний правый
            {-a,  b}  // верхний левый
        };

        // Поворот вершин на угол rect->angle
        for (int i = 0; i < 4; ++i) {
            double local_x = localCorners[i][0];
            double local_y = localCorners[i][1];
            corners[i][0] = x + (local_x * cos(angleRad) - local_y * sin(angleRad)); // поворот вокруг центра
            corners[i][1] = y + (local_x * sin(angleRad) + local_y * cos(angleRad)); // поворот вокруг центра
        }
    }
    bool rectangleCircleIntersect(Rectangle* rect,Circle* circle) {
        double corners[4][2];  // Углы прямоугольника
        getCorners(rect, corners);  // Вычисляем углы прямоугольника

        // Центр окружности
        double cx = circle->geometry[0];
        double cy = circle->geometry[1];
        double radius = circle->geometry[2];

        // Находим минимальные и максимальные координаты прямоугольника по X и Y
        double minX = corners[0][0], maxX = corners[0][0];
        double minY = corners[0][1], maxY = corners[0][1];

        for (int i = 1; i < 4; ++i) {
            if (corners[i][0] < minX) minX = corners[i][0];
            if (corners[i][0] > maxX) maxX = corners[i][0];
            if (corners[i][1] < minY) minY = corners[i][1];
            if (corners[i][1] > maxY) maxY = corners[i][1];
        }

        // Находим ближайшую точку прямоугольника к центру окружности
        double nearestX = cx;
        double nearestY = cy;

        // Проверка по оси X (проекция на стороны прямоугольника)
        if (cx < minX) {
            nearestX = minX;
        }
        else if (cx > maxX) {
            nearestX = maxX;
        }

        // Проверка по оси Y (проекция на стороны прямоугольника)
        if (cy < minY) {
            nearestY = minY;
        }
        else if (cy > maxY) {
            nearestY = maxY;
        }

        // Вычисляем расстояние от ближайшей точки до центра окружности
        double distanceX = cx - nearestX;
        double distanceY = cy - nearestY;
        double distanceSquared = distanceX * distanceX + distanceY * distanceY;

        // Проверка, пересекается ли окружность с прямоугольником
        return distanceSquared <= radius * radius;
    }
    
    bool doProjectionsOverlap(double proj1[2], double proj2[2]) {
        return !(proj1[1] < proj2[0] || proj2[1] < proj1[0]);
    }
    void projectOntoAxis(double axis[2], double corners[4][2], double projection[2]) {
        double min_proj = (corners[0][0] * axis[0] + corners[0][1] * axis[1]);
        double max_proj = min_proj;
        for (int i = 1; i < 4; ++i) {
            double proj = (corners[i][0] * axis[0] + corners[i][1] * axis[1]);
            if (proj < min_proj) min_proj = proj;
            if (proj > max_proj) max_proj = proj;
        }
        projection[0] = min_proj;
        projection[1] = max_proj;
    }
    bool rectanglesIntersect(Rectangle* r1, Rectangle* r2) {
        // Получаем углы для обоих прямоугольников
        double r1_corners[4][2], r2_corners[4][2];
        getCorners(r1, r1_corners);
        getCorners(r2, r2_corners);

        // Перечень осей для проверки
        double axes[4][2] = {
            {r1_corners[1][0] - r1_corners[0][0], r1_corners[1][1] - r1_corners[0][1]}, // ось первой стороны r1
            {r1_corners[3][0] - r1_corners[0][0], r1_corners[3][1] - r1_corners[0][1]}, // ось второй стороны r1
            {r2_corners[1][0] - r2_corners[0][0], r2_corners[1][1] - r2_corners[0][1]}, // ось первой стороны r2
            {r2_corners[3][0] - r2_corners[0][0], r2_corners[3][1] - r2_corners[0][1]}  // ось второй стороны r2
        };

        // Проверяем пересечение проекций на каждой оси
        for (int i = 0; i < 4; ++i) {
            double axis[2] = { axes[i][0], axes[i][1] };
            double r1_projection[2], r2_projection[2];
            projectOntoAxis(axis, r1_corners, r1_projection);
            projectOntoAxis(axis, r2_corners, r2_projection);

            if (!doProjectionsOverlap(r1_projection, r2_projection)) {
                return false; // если проекции не пересекаются, прямоугольники не пересекаются
            }
        }

        return true; // если все проекции пересекаются, прямоугольники пересекаются
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
    Field field = Field(20,0.7);
    field.printFigures();
    field.squreAllFiguresPrint();
    field.combine();
    field.printFigures();
    field.squreAllFiguresPrint();
}


