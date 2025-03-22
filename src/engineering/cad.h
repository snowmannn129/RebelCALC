#ifndef REBELCALC_ENGINEERING_CAD_H
#define REBELCALC_ENGINEERING_CAD_H

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>
#include "../backend/matrix.h"

namespace RebelCalc {
using rebelcalc::Matrix;
namespace Engineering {
namespace CAD {

/**
 * @brief A 3D point/vector class for CAD calculations
 */
class Point3D {
public:
    Point3D() : x(0.0), y(0.0), z(0.0) {}
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
    
    // Basic vector operations
    Point3D operator+(const Point3D& other) const;
    Point3D operator-(const Point3D& other) const;
    Point3D operator*(double scalar) const;
    Point3D operator/(double scalar) const;
    
    // Dot product
    double dot(const Point3D& other) const;
    
    // Cross product
    Point3D cross(const Point3D& other) const;
    
    // Magnitude/length
    double magnitude() const;
    
    // Normalize
    Point3D normalize() const;
    
    // Distance between points
    double distance(const Point3D& other) const;
    
    // Data members
    double x, y, z;
};

/**
 * @brief A line segment in 3D space
 */
class LineSegment {
public:
    LineSegment(const Point3D& start, const Point3D& end) : start(start), end(end) {}
    
    // Length of the line segment
    double length() const;
    
    // Midpoint of the line segment
    Point3D midpoint() const;
    
    // Direction vector (normalized)
    Point3D direction() const;
    
    // Point at parameter t (0 <= t <= 1)
    Point3D pointAt(double t) const;
    
    // Distance from a point to the line segment
    double distanceToPoint(const Point3D& point) const;
    
    // Data members
    Point3D start, end;
};

/**
 * @brief A plane in 3D space defined by a point and normal vector
 */
class Plane {
public:
    Plane(const Point3D& point, const Point3D& normal);
    
    // Check if a point is on the plane
    bool containsPoint(const Point3D& point, double tolerance = 1e-10) const;
    
    // Distance from a point to the plane (signed)
    double distanceToPoint(const Point3D& point) const;
    
    // Project a point onto the plane
    Point3D projectPoint(const Point3D& point) const;
    
    // Intersection of a line with the plane
    // Returns true if there is an intersection, false if line is parallel to plane
    // The intersection point is stored in the 'intersection' parameter
    bool intersectLine(const LineSegment& line, Point3D& intersection) const;
    
    // Data members
    Point3D point;  // A point on the plane
    Point3D normal; // Normal vector to the plane (normalized)
    double d;       // d in the plane equation ax + by + cz + d = 0
};

/**
 * @brief A circle in 3D space
 */
class Circle {
public:
    Circle(const Point3D& center, const Point3D& normal, double radius);
    
    // Area of the circle
    double area() const;
    
    // Circumference of the circle
    double circumference() const;
    
    // Check if a point is on the circle
    bool containsPoint(const Point3D& point, double tolerance = 1e-10) const;
    
    // Point on the circle at angle theta (in radians)
    Point3D pointAtAngle(double theta) const;
    
    // Data members
    Point3D center;  // Center of the circle
    Point3D normal;  // Normal vector to the plane containing the circle (normalized)
    double radius;   // Radius of the circle
    
private:
    // Cached basis vectors for the circle's plane
    Point3D u, v;
    
    // Initialize the basis vectors
    void initBasis();
};

/**
 * @brief A triangle in 3D space
 */
class Triangle {
public:
    Triangle(const Point3D& p1, const Point3D& p2, const Point3D& p3);
    
    // Area of the triangle
    double area() const;
    
    // Perimeter of the triangle
    double perimeter() const;
    
    // Normal vector to the triangle
    Point3D normal() const;
    
    // Check if a point is inside the triangle
    bool containsPoint(const Point3D& point, double tolerance = 1e-10) const;
    
    // Barycentric coordinates of a point with respect to the triangle
    std::array<double, 3> barycentricCoordinates(const Point3D& point) const;
    
    // Data members
    Point3D p1, p2, p3;
};

/**
 * @brief A polygon in 3D space
 */
class Polygon {
public:
    Polygon(const std::vector<Point3D>& vertices);
    
    // Area of the polygon (assumes it's planar and simple)
    double area() const;
    
    // Perimeter of the polygon
    double perimeter() const;
    
    // Check if the polygon is convex
    bool isConvex() const;
    
    // Check if a point is inside the polygon (assumes it's planar)
    bool containsPoint(const Point3D& point, double tolerance = 1e-10) const;
    
    // Triangulate the polygon
    std::vector<Triangle> triangulate() const;
    
    // Data members
    std::vector<Point3D> vertices;
};

/**
 * @brief Calculates the volume of a tetrahedron defined by four points
 * 
 * @param p1 First point
 * @param p2 Second point
 * @param p3 Third point
 * @param p4 Fourth point
 * @return Volume of the tetrahedron
 */
double tetrahedronVolume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4);

/**
 * @brief Calculates the centroid of a set of points
 * 
 * @param points Vector of points
 * @return Centroid point
 */
Point3D centroid(const std::vector<Point3D>& points);

/**
 * @brief Calculates the best-fit plane for a set of points
 * 
 * @param points Vector of points
 * @return Best-fit plane
 */
Plane bestFitPlane(const std::vector<Point3D>& points);

/**
 * @brief Calculates the convex hull of a set of points in 3D (using the Gift Wrapping algorithm)
 * 
 * @param points Vector of points
 * @return Vector of triangles forming the convex hull
 */
std::vector<Triangle> convexHull3D(const std::vector<Point3D>& points);

/**
 * @brief Calculates the intersection of two planes
 * 
 * @param plane1 First plane
 * @param plane2 Second plane
 * @param line Line of intersection (output parameter)
 * @return true if the planes intersect, false if they are parallel
 */
bool planePlaneIntersection(const Plane& plane1, const Plane& plane2, LineSegment& line);

/**
 * @brief Calculates the intersection of three planes
 * 
 * @param plane1 First plane
 * @param plane2 Second plane
 * @param plane3 Third plane
 * @param point Point of intersection (output parameter)
 * @return true if the planes intersect at a single point, false otherwise
 */
bool threePlaneIntersection(const Plane& plane1, const Plane& plane2, const Plane& plane3, Point3D& point);

/**
 * @brief Calculates the minimum distance between two line segments in 3D
 * 
 * @param line1 First line segment
 * @param line2 Second line segment
 * @return Minimum distance between the line segments
 */
double lineLineDistance(const LineSegment& line1, const LineSegment& line2);

/**
 * @brief Calculates the minimum distance between a point and a triangle
 * 
 * @param point Point
 * @param triangle Triangle
 * @return Minimum distance between the point and the triangle
 */
double pointTriangleDistance(const Point3D& point, const Triangle& triangle);

/**
 * @brief Calculates the intersection of a line and a triangle
 * 
 * @param line Line segment
 * @param triangle Triangle
 * @param intersection Point of intersection (output parameter)
 * @return true if the line intersects the triangle, false otherwise
 */
bool lineTriangleIntersection(const LineSegment& line, const Triangle& triangle, Point3D& intersection);

/**
 * @brief Calculates the intersection of a ray and a triangle (Möller–Trumbore algorithm)
 * 
 * @param origin Ray origin
 * @param direction Ray direction (normalized)
 * @param triangle Triangle
 * @param t Distance along the ray to the intersection (output parameter)
 * @param intersection Point of intersection (output parameter)
 * @return true if the ray intersects the triangle, false otherwise
 */
bool rayTriangleIntersection(const Point3D& origin, const Point3D& direction, 
                             const Triangle& triangle, double& t, Point3D& intersection);

/**
 * @brief Calculates the bounding box of a set of points
 * 
 * @param points Vector of points
 * @return std::array containing the minimum and maximum points of the bounding box
 */
std::array<Point3D, 2> boundingBox(const std::vector<Point3D>& points);

/**
 * @brief Calculates the bounding sphere of a set of points
 * 
 * @param points Vector of points
 * @param center Center of the bounding sphere (output parameter)
 * @param radius Radius of the bounding sphere (output parameter)
 */
void boundingSphere(const std::vector<Point3D>& points, Point3D& center, double& radius);

/**
 * @brief Calculates the moment of inertia tensor for a set of points with equal mass
 * 
 * @param points Vector of points
 * @param mass Total mass of the system
 * @return 3x3 inertia tensor matrix
 */
rebelcalc::Matrix inertiaTensor(const std::vector<Point3D>& points, double mass);

/**
 * @brief Calculates the principal axes and moments of inertia
 * 
 * @param inertiaTensor 3x3 inertia tensor matrix
 * @param principalAxes 3x3 matrix where columns are the principal axes (output parameter)
 * @param principalMoments Vector of principal moments of inertia (output parameter)
 */
void principalAxesAndMoments(const rebelcalc::Matrix& inertiaTensor, rebelcalc::Matrix& principalAxes, std::vector<double>& principalMoments);

} // namespace CAD
} // namespace Engineering
} // namespace RebelCalc

#endif // REBELCALC_ENGINEERING_CAD_H
