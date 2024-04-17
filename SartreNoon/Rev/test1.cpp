#include <iostream>
#include <vector>

using namespace std;

struct Point {
    int x;
    int y;
};

int main() {
    // Using custom struct
    vector<Point> points;

    // Adding points to the vector
    points.push_back({1, 2});
    points.push_back({3, 4});
    points.push_back({5, 6});

    // Accessing points
    for (const auto& point : points) {
        cout << "X: " << point.x << ", Y: " << point.y << endl;
    }

    return 0;
}
