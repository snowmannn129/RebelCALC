#include "cad.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <cmath>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RebelCalc {
namespace Engineering {
namespace CAD {

//
// Point3D implementation
//

Point3D Point3D::operator+(const Point3D& other) const {
    return Point3D(x + other.x, y + other.y, z + other.z);
}

Point3D Point3D::operator-(const Point3D& other) const {
    return Point3D(x - other.x, y - other.y, z - other.z);
}

Point3D Point3D::operator*(double scalar) const {
    return Point3D(x * scalar, y * scalar, z * scalar);
}

Point3D Point3D::operator/(double scalar) const {
    if (std::abs(scalar) < 1e-10) {
        throw std::invalid_argument("Division by zero or near-zero value");
    }
    return Point3D(x / scalar, y / scalar, z / scalar);
}

double Point3D::dot(const Point3D& other) const {
    return x * other.x + y * other.y + z * other.z;
}

Point3D Point3D::cross(const Point3D& other) const {
    return Point3D(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

double Point3D::magnitude() const {
    return std::sqrt(x * x + y * y + z * z);
}

Point3D Point3D::normalize() const {
    double mag = magnitude();
    if (mag < 1e-10) {
        throw std::runtime_error("Cannot normalize a zero vector");
    }
    return *this / mag;
}

double Point3D::distance(const Point3D& other) const {
    return (*this - other).magnitude();
}

//
// LineSegment implementation
//

double LineSegment::length() const {
    return start.distance(end);
}

Point3D LineSegment::midpoint() const {
    return Point3D(
        (start.x + end.x) / 2.0,
        (start.y + end.y) / 2.0,
        (start.z + end.z) / 2.0
    );
}

Point3D LineSegment::direction() const {
    return (end - start).normalize();
}

Point3D LineSegment::pointAt(double t) const {
    if (t < 0.0 || t > 1.0) {
        throw std::out_of_range("Parameter t must be between 0 and 1");
    }
    return start + (end - start) * t;
}

double LineSegment::distanceToPoint(const Point3D& point) const {
    Point3D v = end - start;
    Point3D w = point - start;
    
    double c1 = w.dot(v);
    if (c1 <= 0) {
        return w.magnitude();
    }
    
    double c2 = v.dot(v);
    if (c2 <= c1) {
        return (point - end).magnitude();
    }
    
    double b = c1 / c2;
    Point3D pb = start + v * b;
    return (point - pb).magnitude();
}

//
// Plane implementation
//

Plane::Plane(const Point3D& point, const Point3D& normal) : point(point) {
    this->normal = normal.normalize();
    d = -this->normal.dot(point);
}

bool Plane::containsPoint(const Point3D& point, double tolerance) const {
    return std::abs(distanceToPoint(point)) < tolerance;
}

double Plane::distanceToPoint(const Point3D& point) const {
    return normal.dot(point) + d;
}

Point3D Plane::projectPoint(const Point3D& point) const {
    double dist = distanceToPoint(point);
    return point - normal * dist;
}

bool Plane::intersectLine(const LineSegment& line, Point3D& intersection) const {
    Point3D dir = line.end - line.start;
    double denom = normal.dot(dir);
    
    // Check if line is parallel to plane
    if (std::abs(denom) < 1e-10) {
        return false;
    }
    
    double t = -(normal.dot(line.start) + d) / denom;
    
    // Check if intersection is within the line segment
    if (t < 0.0 || t > 1.0) {
        return false;
    }
    
    intersection = line.start + dir * t;
    return true;
}

//
// Circle implementation
//

Circle::Circle(const Point3D& center, const Point3D& normal, double radius)
    : center(center), radius(radius) {
    this->normal = normal.normalize();
    initBasis();
}

void Circle::initBasis() {
    // Create two orthogonal vectors in the plane of the circle
    Point3D arbitrary;
    if (std::abs(this->normal.x) < std::abs(this->normal.y) && 
        std::abs(this->normal.x) < std::abs(this->normal.z)) {
        arbitrary = Point3D(1, 0, 0);
    } else if (std::abs(this->normal.y) < std::abs(this->normal.z)) {
        arbitrary = Point3D(0, 1, 0);
    } else {
        arbitrary = Point3D(0, 0, 1);
    }
    
    u = this->normal.cross(arbitrary).normalize();
    v = this->normal.cross(u).normalize();
}

double Circle::area() const {
    return M_PI * radius * radius;
}

double Circle::circumference() const {
    return 2.0 * M_PI * radius;
}

bool Circle::containsPoint(const Point3D& point, double tolerance) const {
    // Check if point is on the plane of the circle
    Plane plane(center, normal);
    if (!plane.containsPoint(point, tolerance)) {
        return false;
    }
    
    // Check if point is at the right distance from the center
    return std::abs(center.distance(point) - radius) < tolerance;
}

Point3D Circle::pointAtAngle(double theta) const {
    return center + u * (radius * std::cos(theta)) + v * (radius * std::sin(theta));
}

//
// Triangle implementation
//

Triangle::Triangle(const Point3D& p1, const Point3D& p2, const Point3D& p3)
    : p1(p1), p2(p2), p3(p3) {}

double Triangle::area() const {
    Point3D v1 = p2 - p1;
    Point3D v2 = p3 - p1;
    return 0.5 * v1.cross(v2).magnitude();
}

double Triangle::perimeter() const {
    return p1.distance(p2) + p2.distance(p3) + p3.distance(p1);
}

Point3D Triangle::normal() const {
    Point3D v1 = p2 - p1;
    Point3D v2 = p3 - p1;
    return v1.cross(v2).normalize();
}

bool Triangle::containsPoint(const Point3D& point, double tolerance) const {
    // Check if point is on the plane of the triangle
    Plane plane(p1, normal());
    if (!plane.containsPoint(point, tolerance)) {
        return false;
    }
    
    // Use barycentric coordinates to check if point is inside the triangle
    std::array<double, 3> bary = barycentricCoordinates(point);
    return bary[0] >= -tolerance && bary[1] >= -tolerance && bary[2] >= -tolerance &&
           bary[0] <= 1.0 + tolerance && bary[1] <= 1.0 + tolerance && bary[2] <= 1.0 + tolerance;
}

std::array<double, 3> Triangle::barycentricCoordinates(const Point3D& point) const {
    Point3D v0 = p2 - p1;
    Point3D v1 = p3 - p1;
    Point3D v2 = point - p1;
    
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;
    
    return {u, v, w};
}

//
// Polygon implementation
//

Polygon::Polygon(const std::vector<Point3D>& vertices) : vertices(vertices) {
    if (vertices.size() < 3) {
        throw std::invalid_argument("A polygon must have at least 3 vertices");
    }
}

double Polygon::area() const {
    // Compute area using the shoelace formula (for planar polygons)
    // First, project the polygon onto a plane
    Point3D normal = bestFitPlane(vertices).normal;
    
    // Find the axis with the largest normal component
    int axis = 0;
    if (std::abs(normal.y) > std::abs(normal.x) && std::abs(normal.y) > std::abs(normal.z)) {
        axis = 1;
    } else if (std::abs(normal.z) > std::abs(normal.x) && std::abs(normal.z) > std::abs(normal.y)) {
        axis = 2;
    }
    
    // Project onto the plane perpendicular to the chosen axis
    double area = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Point3D& v1 = vertices[i];
        const Point3D& v2 = vertices[(i + 1) % vertices.size()];
        
        if (axis == 0) {
            area += v1.y * v2.z - v1.z * v2.y;
        } else if (axis == 1) {
            area += v1.z * v2.x - v1.x * v2.z;
        } else {
            area += v1.x * v2.y - v1.y * v2.x;
        }
    }
    
    return 0.5 * std::abs(area);
}

double Polygon::perimeter() const {
    double perim = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        perim += vertices[i].distance(vertices[(i + 1) % vertices.size()]);
    }
    return perim;
}

bool Polygon::isConvex() const {
    if (vertices.size() < 3) {
        return false;
    }
    
    // For a 3D polygon, we need to ensure all vertices lie in the same plane
    Plane plane = bestFitPlane(vertices);
    
    // Check if all vertices lie on the plane
    for (const auto& vertex : vertices) {
        if (!plane.containsPoint(vertex)) {
            return false;
        }
    }
    
    // Project the polygon onto the plane
    std::vector<Point3D> projectedVertices;
    for (const auto& vertex : vertices) {
        projectedVertices.push_back(plane.projectPoint(vertex));
    }
    
    // Find the axis with the largest normal component
    int axis = 0;
    if (std::abs(plane.normal.y) > std::abs(plane.normal.x) && 
        std::abs(plane.normal.y) > std::abs(plane.normal.z)) {
        axis = 1;
    } else if (std::abs(plane.normal.z) > std::abs(plane.normal.x) && 
               std::abs(plane.normal.z) > std::abs(plane.normal.y)) {
        axis = 2;
    }
    
    // Check if the polygon is convex by checking the sign of the cross product
    bool sign = false;
    bool signSet = false;
    
    for (size_t i = 0; i < projectedVertices.size(); ++i) {
        const Point3D& v1 = projectedVertices[i];
        const Point3D& v2 = projectedVertices[(i + 1) % projectedVertices.size()];
        const Point3D& v3 = projectedVertices[(i + 2) % projectedVertices.size()];
        
        Point3D edge1, edge2;
        if (axis == 0) {
            edge1 = Point3D(0, v2.y - v1.y, v2.z - v1.z);
            edge2 = Point3D(0, v3.y - v2.y, v3.z - v2.z);
        } else if (axis == 1) {
            edge1 = Point3D(v2.x - v1.x, 0, v2.z - v1.z);
            edge2 = Point3D(v3.x - v2.x, 0, v3.z - v2.z);
        } else {
            edge1 = Point3D(v2.x - v1.x, v2.y - v1.y, 0);
            edge2 = Point3D(v3.x - v2.x, v3.y - v2.y, 0);
        }
        
        Point3D cross = edge1.cross(edge2);
        double crossProduct;
        if (axis == 0) {
            crossProduct = cross.x;
        } else if (axis == 1) {
            crossProduct = cross.y;
        } else {
            crossProduct = cross.z;
        }
        
        bool currentSign = crossProduct > 0;
        if (!signSet) {
            sign = currentSign;
            signSet = true;
        } else if (sign != currentSign) {
            return false;
        }
    }
    
    return true;
}

bool Polygon::containsPoint(const Point3D& point, double tolerance) const {
    // Check if point is on the plane of the polygon
    Plane plane = bestFitPlane(vertices);
    if (!plane.containsPoint(point, tolerance)) {
        return false;
    }
    
    // Project the polygon and the point onto a 2D plane
    int axis = 0;
    if (std::abs(plane.normal.y) > std::abs(plane.normal.x) && 
        std::abs(plane.normal.y) > std::abs(plane.normal.z)) {
        axis = 1;
    } else if (std::abs(plane.normal.z) > std::abs(plane.normal.x) && 
               std::abs(plane.normal.z) > std::abs(plane.normal.y)) {
        axis = 2;
    }
    
    // Use the ray casting algorithm to check if the point is inside the polygon
    bool inside = false;
    for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
        double xi, yi, xj, yj, xp, yp;
        
        if (axis == 0) {
            xi = vertices[i].y; yi = vertices[i].z;
            xj = vertices[j].y; yj = vertices[j].z;
            xp = point.y; yp = point.z;
        } else if (axis == 1) {
            xi = vertices[i].x; yi = vertices[i].z;
            xj = vertices[j].x; yj = vertices[j].z;
            xp = point.x; yp = point.z;
        } else {
            xi = vertices[i].x; yi = vertices[i].y;
            xj = vertices[j].x; yj = vertices[j].y;
            xp = point.x; yp = point.y;
        }
        
        if (((yi > yp) != (yj > yp)) &&
            (xp < (xj - xi) * (yp - yi) / (yj - yi) + xi)) {
            inside = !inside;
        }
    }
    
    return inside;
}

std::vector<Triangle> Polygon::triangulate() const {
    // Simple ear clipping triangulation for convex polygons
    if (vertices.size() < 3) {
        return {};
    }
    
    if (vertices.size() == 3) {
        return {Triangle(vertices[0], vertices[1], vertices[2])};
    }
    
    std::vector<Triangle> triangles;
    std::vector<Point3D> remainingVertices = vertices;
    
    while (remainingVertices.size() > 3) {
        bool earFound = false;
        
        for (size_t i = 0; i < remainingVertices.size(); ++i) {
            size_t prev = (i + remainingVertices.size() - 1) % remainingVertices.size();
            size_t next = (i + 1) % remainingVertices.size();
            
            Triangle ear(remainingVertices[prev], remainingVertices[i], remainingVertices[next]);
            
            // Check if this is a valid ear (no other vertices inside)
            bool isEar = true;
            for (size_t j = 0; j < remainingVertices.size(); ++j) {
                if (j != prev && j != i && j != next && 
                    ear.containsPoint(remainingVertices[j])) {
                    isEar = false;
                    break;
                }
            }
            
            if (isEar) {
                triangles.push_back(ear);
                remainingVertices.erase(remainingVertices.begin() + i);
                earFound = true;
                break;
            }
        }
        
        if (!earFound) {
            // If no ear is found, the polygon might be non-simple or have collinear points
            // In this case, just create a triangle with consecutive vertices
            triangles.push_back(Triangle(remainingVertices[0], remainingVertices[1], remainingVertices[2]));
            remainingVertices.erase(remainingVertices.begin() + 1);
        }
    }
    
    // Add the last triangle
    triangles.push_back(Triangle(remainingVertices[0], remainingVertices[1], remainingVertices[2]));
    
    return triangles;
}

//
// Utility functions
//

double tetrahedronVolume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4) {
    Point3D v1 = p2 - p1;
    Point3D v2 = p3 - p1;
    Point3D v3 = p4 - p1;
    
    // Volume = (1/6) * |triple scalar product|
    return std::abs(v1.dot(v2.cross(v3))) / 6.0;
}

Point3D centroid(const std::vector<Point3D>& points) {
    if (points.empty()) {
        throw std::invalid_argument("Cannot compute centroid of empty point set");
    }
    
    Point3D sum(0, 0, 0);
    for (const auto& point : points) {
        sum = sum + point;
    }
    
    return sum / static_cast<double>(points.size());
}

Plane bestFitPlane(const std::vector<Point3D>& points) {
    if (points.size() < 3) {
        throw std::invalid_argument("At least 3 points are required to fit a plane");
    }
    
    // Compute the centroid
    Point3D c = centroid(points);
    
    // Compute the covariance matrix
    double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    
    for (const auto& p : points) {
        Point3D pc = p - c;
        xx += pc.x * pc.x;
        xy += pc.x * pc.y;
        xz += pc.x * pc.z;
        yy += pc.y * pc.y;
        yz += pc.y * pc.z;
        zz += pc.z * pc.z;
    }
    
    // Find the smallest eigenvalue and corresponding eigenvector
    // This is a simplified approach for a 3x3 matrix
    rebelcalc::Matrix cov(3, 3);
    cov(0, 0) = xx; cov(0, 1) = xy; cov(0, 2) = xz;
    cov(1, 0) = xy; cov(1, 1) = yy; cov(1, 2) = yz;
    cov(2, 0) = xz; cov(2, 1) = yz; cov(2, 2) = zz;
    
    // Use power iteration to find the eigenvector with the smallest eigenvalue
    // (This is a simplified approach and might not be the most accurate)
    Point3D normal(1, 1, 1);
    for (int i = 0; i < 10; ++i) {
        double nx = cov(0, 0) * normal.x + cov(0, 1) * normal.y + cov(0, 2) * normal.z;
        double ny = cov(1, 0) * normal.x + cov(1, 1) * normal.y + cov(1, 2) * normal.z;
        double nz = cov(2, 0) * normal.x + cov(2, 1) * normal.y + cov(2, 2) * normal.z;
        
        normal = Point3D(nx, ny, nz).normalize();
    }
    
    return Plane(c, normal);
}

std::vector<Triangle> convexHull3D(const std::vector<Point3D>& points) {
    if (points.size() < 4) {
        throw std::invalid_argument("At least 4 non-coplanar points are required for a 3D convex hull");
    }
    
    // This is a simplified implementation of the Gift Wrapping algorithm
    // for 3D convex hull. A more robust implementation would use a more
    // sophisticated algorithm like QuickHull or Incremental Hull.
    
    // Find the point with the minimum x-coordinate
    auto minXIt = std::min_element(points.begin(), points.end(),
                                  [](const Point3D& a, const Point3D& b) {
                                      return a.x < b.x;
                                  });
    Point3D minX = *minXIt;
    
    // Find the point farthest from minX
    auto maxDistIt = std::max_element(points.begin(), points.end(),
                                     [&minX](const Point3D& a, const Point3D& b) {
                                         return minX.distance(a) < minX.distance(b);
                                     });
    Point3D maxDist = *maxDistIt;
    
    // Find the point that forms the largest area with minX and maxDist
    auto maxAreaIt = std::max_element(points.begin(), points.end(),
                                     [&minX, &maxDist](const Point3D& a, const Point3D& b) {
                                         Point3D v1 = maxDist - minX;
                                         Point3D v2a = a - minX;
                                         Point3D v2b = b - minX;
                                         return v1.cross(v2a).magnitude() < v1.cross(v2b).magnitude();
                                     });
    Point3D maxArea = *maxAreaIt;
    
    // Initial triangle
    Triangle initialTriangle(minX, maxDist, maxArea);
    
    // Find the point with the maximum distance from the initial triangle
    auto maxDistFromTriangleIt = std::max_element(points.begin(), points.end(),
                                                 [&initialTriangle](const Point3D& a, const Point3D& b) {
                                                     Plane plane(initialTriangle.p1, initialTriangle.normal());
                                                     return std::abs(plane.distanceToPoint(a)) < 
                                                            std::abs(plane.distanceToPoint(b));
                                                 });
    Point3D apex = *maxDistFromTriangleIt;
    
    // Initial tetrahedron
    std::vector<Triangle> hull = {
        Triangle(minX, maxDist, maxArea),
        Triangle(minX, maxDist, apex),
        Triangle(minX, maxArea, apex),
        Triangle(maxDist, maxArea, apex)
    };
    
    // Ensure all triangles have outward-facing normals
    for (auto& triangle : hull) {
        Point3D centroid = (triangle.p1 + triangle.p2 + triangle.p3) / 3.0;
        Point3D normal = triangle.normal();
        Point3D outward = centroid - apex;
        if (normal.dot(outward) < 0) {
            // Swap two vertices to flip the normal
            std::swap(triangle.p2, triangle.p3);
        }
    }
    
    // This is a simplified implementation and doesn't handle all cases
    // A complete implementation would iteratively add points to the hull
    
    return hull;
}

bool planePlaneIntersection(const Plane& plane1, const Plane& plane2, LineSegment& line) {
    // Direction of the intersection line is the cross product of the normals
    Point3D direction = plane1.normal.cross(plane2.normal);
    
    // Check if planes are parallel
    if (direction.magnitude() < 1e-10) {
        return false;
    }
    
    // Normalize the direction
    direction = direction.normalize();
    
    // Find a point on the intersection line
    // We solve the system of equations:
    // n1 . p + d1 = 0
    // n2 . p + d2 = 0
    // where n1, n2 are the normals, p is a point on the line, and d1, d2 are the plane constants
    
    // Find the component with the largest absolute value in the direction
    int maxComponent = 0;
    if (std::abs(direction.y) > std::abs(direction.x) && 
        std::abs(direction.y) > std::abs(direction.z)) {
        maxComponent = 1;
    } else if (std::abs(direction.z) > std::abs(direction.x) && 
               std::abs(direction.z) > std::abs(direction.y)) {
        maxComponent = 2;
    }
    
    // Set that component to 0 and solve for the other two
    Point3D point;
    if (maxComponent == 0) {
        point.x = 0;
        // Solve for y and z
        double a1 = plane1.normal.y, b1 = plane1.normal.z, c1 = -plane1.d;
        double a2 = plane2.normal.y, b2 = plane2.normal.z, c2 = -plane2.d;
        
        double det = a1 * b2 - a2 * b1;
        if (std::abs(det) < 1e-10) {
            // The system is singular, try another approach
            point.y = 0;
            point.z = -plane1.d / plane1.normal.z;
        } else {
            point.y = (c1 * b2 - c2 * b1) / det;
            point.z = (a1 * c2 - a2 * c1) / det;
        }
    } else if (maxComponent == 1) {
        point.y = 0;
        // Solve for x and z
        double a1 = plane1.normal.x, b1 = plane1.normal.z, c1 = -plane1.d;
        double a2 = plane2.normal.x, b2 = plane2.normal.z, c2 = -plane2.d;
        
        double det = a1 * b2 - a2 * b1;
        if (std::abs(det) < 1e-10) {
            // The system is singular, try another approach
            point.x = 0;
            point.z = -plane1.d / plane1.normal.z;
        } else {
            point.x = (c1 * b2 - c2 * b1) / det;
            point.z = (a1 * c2 - a2 * c1) / det;
        }
    } else {
        point.z = 0;
        // Solve for x and y
        double a1 = plane1.normal.x, b1 = plane1.normal.y, c1 = -plane1.d;
        double a2 = plane2.normal.x, b2 = plane2.normal.y, c2 = -plane2.d;
        
        double det = a1 * b2 - a2 * b1;
        if (std::abs(det) < 1e-10) {
            // The system is singular, try another approach
            point.x = 0;
            point.y = -plane1.d / plane1.normal.y;
        } else {
            point.x = (c1 * b2 - c2 * b1) / det;
            point.y = (a1 * c2 - a2 * c1) / det;
        }
    }
    
    // Create a line segment along the intersection line
    // We use an arbitrary length for the line segment
    double length = 10.0;
    line = LineSegment(point - direction * length, point + direction * length);
    
    return true;
}

bool threePlaneIntersection(const Plane& plane1, const Plane& plane2, const Plane& plane3, Point3D& point) {
    // Create the system of equations:
    // n1 . p + d1 = 0
    // n2 . p + d2 = 0
    // n3 . p + d3 = 0
    // where n1, n2, n3 are the normals, p is the intersection point, and d1, d2, d3 are the plane constants
    
    // Create the coefficient matrix
    rebelcalc::Matrix A(3, 3);
    A(0, 0) = plane1.normal.x; A(0, 1) = plane1.normal.y; A(0, 2) = plane1.normal.z;
    A(1, 0) = plane2.normal.x; A(1, 1) = plane2.normal.y; A(1, 2) = plane2.normal.z;
    A(2, 0) = plane3.normal.x; A(2, 1) = plane3.normal.y; A(2, 2) = plane3.normal.z;
    
    // Create the right-hand side vector
    std::vector<double> b = {-plane1.d, -plane2.d, -plane3.d};
    
    // Check if the system has a unique solution
    double det = A.determinant();
    if (std::abs(det) < 1e-10) {
        return false;
    }
    
    // Solve the system
    std::vector<double> solution = A.solve(b);
    point = Point3D(solution[0], solution[1], solution[2]);
    
    return true;
}

double lineLineDistance(const LineSegment& line1, const LineSegment& line2) {
    Point3D u = line1.end - line1.start;
    Point3D v = line2.end - line2.start;
    Point3D w = line1.start - line2.start;
    
    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    
    double D = a * c - b * b;
    
    // Lines are parallel
    if (D < 1e-10) {
        // Find the distance between line1.start and line2
        return line2.distanceToPoint(line1.start);
    }
    
    double sc = (b * e - c * d) / D;
    double tc = (a * e - b * d) / D;
    
    // Clamp sc and tc to [0, 1]
    sc = std::max(0.0, std::min(1.0, sc));
    tc = std::max(0.0, std::min(1.0, tc));
    
    // Compute the closest points
    Point3D p1 = line1.start + u * sc;
    Point3D p2 = line2.start + v * tc;
    
    return p1.distance(p2);
}

double pointTriangleDistance(const Point3D& point, const Triangle& triangle) {
    // Project the point onto the plane of the triangle
    Plane plane(triangle.p1, triangle.normal());
    Point3D projected = plane.projectPoint(point);
    
    // Check if the projected point is inside the triangle
    if (triangle.containsPoint(projected)) {
        // Distance is just the distance to the plane
        return std::abs(plane.distanceToPoint(point));
    }
    
    // Find the closest point on the edges of the triangle
    LineSegment edge1(triangle.p1, triangle.p2);
    LineSegment edge2(triangle.p2, triangle.p3);
    LineSegment edge3(triangle.p3, triangle.p1);
    
    double d1 = edge1.distanceToPoint(point);
    double d2 = edge2.distanceToPoint(point);
    double d3 = edge3.distanceToPoint(point);
    
    return std::min({d1, d2, d3});
}

bool lineTriangleIntersection(const LineSegment& line, const Triangle& triangle, Point3D& intersection) {
    // Check if the line intersects the plane of the triangle
    Plane plane(triangle.p1, triangle.normal());
    if (!plane.intersectLine(line, intersection)) {
        return false;
    }
    
    // Check if the intersection point is inside the triangle
    return triangle.containsPoint(intersection);
}

bool rayTriangleIntersection(const Point3D& origin, const Point3D& direction, 
                             const Triangle& triangle, double& t, Point3D& intersection) {
    // Möller–Trumbore algorithm
    Point3D edge1 = triangle.p2 - triangle.p1;
    Point3D edge2 = triangle.p3 - triangle.p1;
    Point3D h = direction.cross(edge2);
    double a = edge1.dot(h);
    
    // Check if ray is parallel to triangle
    if (std::abs(a) < 1e-10) {
        return false;
    }
    
    double f = 1.0 / a;
    Point3D s = origin - triangle.p1;
    double u = f * s.dot(h);
    
    // Check if intersection is outside the triangle
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    
    Point3D q = s.cross(edge1);
    double v = f * direction.dot(q);
    
    // Check if intersection is outside the triangle
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    
    // Compute the distance along the ray
    t = f * edge2.dot(q);
    
    // Check if intersection is behind the ray origin
    if (t < 0.0) {
        return false;
    }
    
    // Compute the intersection point
    intersection = origin + direction * t;
    
    return true;
}

std::array<Point3D, 2> boundingBox(const std::vector<Point3D>& points) {
    if (points.empty()) {
        throw std::invalid_argument("Cannot compute bounding box of empty point set");
    }
    
    Point3D min = points[0];
    Point3D max = points[0];
    
    for (const auto& point : points) {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);
        
        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
    
    return {min, max};
}

void boundingSphere(const std::vector<Point3D>& points, Point3D& center, double& radius) {
    if (points.empty()) {
        throw std::invalid_argument("Cannot compute bounding sphere of empty point set");
    }
    
    // Ritter's algorithm for approximate bounding sphere
    
    // Start with any point
    center = points[0];
    radius = 0.0;
    
    // Find the point farthest from the current center
    for (int i = 0; i < 2; ++i) {
        Point3D farthest = points[0];
        double maxDist = center.distance(farthest);
        
        for (const auto& point : points) {
            double dist = center.distance(point);
            if (dist > maxDist) {
                maxDist = dist;
                farthest = point;
            }
        }
        
        // Update center and radius
        center = (center + farthest) / 2.0;
        radius = center.distance(farthest);
    }
    
    // Ensure all points are inside the sphere
    for (const auto& point : points) {
        double dist = center.distance(point);
        if (dist > radius) {
            // Expand the sphere to include this point
            double newRadius = (radius + dist) / 2.0;
            double ratio = (newRadius - radius) / dist;
            center = center + (point - center) * ratio;
            radius = newRadius;
        }
    }
}

rebelcalc::Matrix inertiaTensor(const std::vector<Point3D>& points, double mass) {
    if (points.empty()) {
        throw std::invalid_argument("Cannot compute inertia tensor of empty point set");
    }
    
    // Compute the centroid
    Point3D c = centroid(points);
    
    // Compute the inertia tensor
    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
    
    double pointMass = mass / points.size();
    
    for (const auto& p : points) {
        Point3D r = p - c;
        
        Ixx += pointMass * (r.y * r.y + r.z * r.z);
        Iyy += pointMass * (r.x * r.x + r.z * r.z);
        Izz += pointMass * (r.x * r.x + r.y * r.y);
        
        Ixy -= pointMass * r.x * r.y;
        Ixz -= pointMass * r.x * r.z;
        Iyz -= pointMass * r.y * r.z;
    }
    
    rebelcalc::Matrix I(3, 3);
    I(0, 0) = Ixx; I(0, 1) = Ixy; I(0, 2) = Ixz;
    I(1, 0) = Ixy; I(1, 1) = Iyy; I(1, 2) = Iyz;
    I(2, 0) = Ixz; I(2, 1) = Iyz; I(2, 2) = Izz;
    
    return I;
}

void principalAxesAndMoments(const rebelcalc::Matrix& inertiaTensor, rebelcalc::Matrix& principalAxes, std::vector<double>& principalMoments) {
    // This is a simplified approach that assumes the inertia tensor is symmetric
    // A more robust implementation would use a proper eigenvalue decomposition algorithm
    
    // For now, we'll just return the inertia tensor as is
    principalAxes = rebelcalc::Matrix::identity(3);
    principalMoments = {inertiaTensor(0, 0), inertiaTensor(1, 1), inertiaTensor(2, 2)};
}
